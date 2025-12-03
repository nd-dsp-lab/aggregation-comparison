#ifdef __gramine__
#include <gramine/api.hh>
#endif

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <cstring>
#include <openssl/rsa.h>
#include <openssl/pem.h>
#include <openssl/err.h>
#include <openssl/evp.h>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <fstream>

struct EncryptedData {
    uint32_t device_id;
    std::vector<uint8_t> encrypted_value;
};

struct SharedMemory {
    std::vector<EncryptedData> all_data;
    std::vector<uint64_t> round_sums;
    EVP_PKEY* private_key;
    EVP_PKEY* public_key;
};

EVP_PKEY* generate_rsa_keypair() {
    EVP_PKEY_CTX* ctx = EVP_PKEY_CTX_new_id(EVP_PKEY_RSA, NULL);
    if (!ctx) throw std::runtime_error("Failed to create EVP_PKEY_CTX");

    if (EVP_PKEY_keygen_init(ctx) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw std::runtime_error("Failed to init keygen");
    }

    if (EVP_PKEY_CTX_set_rsa_keygen_bits(ctx, 3072) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw std::runtime_error("Failed to set key bits");
    }

    EVP_PKEY* pkey = NULL;
    if (EVP_PKEY_keygen(ctx, &pkey) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw std::runtime_error("Failed to generate keypair");
    }

    EVP_PKEY_CTX_free(ctx);
    return pkey;
}

std::vector<uint8_t> rsa_encrypt(EVP_PKEY* public_key, const uint8_t* plaintext, size_t plaintext_len) {
    EVP_PKEY_CTX* ctx = EVP_PKEY_CTX_new(public_key, NULL);
    if (!ctx) throw std::runtime_error("Failed to create encrypt context");

    if (EVP_PKEY_encrypt_init(ctx) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw std::runtime_error("Failed to init encryption");
    }

    if (EVP_PKEY_CTX_set_rsa_padding(ctx, RSA_PKCS1_OAEP_PADDING) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw std::runtime_error("Failed to set padding");
    }

    size_t outlen;
    if (EVP_PKEY_encrypt(ctx, NULL, &outlen, plaintext, plaintext_len) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw std::runtime_error("Failed to determine ciphertext length");
    }

    std::vector<uint8_t> ciphertext(outlen);
    if (EVP_PKEY_encrypt(ctx, ciphertext.data(), &outlen, plaintext, plaintext_len) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw std::runtime_error("Failed to encrypt");
    }

    ciphertext.resize(outlen);
    EVP_PKEY_CTX_free(ctx);
    return ciphertext;
}

uint64_t rsa_decrypt(EVP_PKEY* private_key, const std::vector<uint8_t>& ciphertext) {
    EVP_PKEY_CTX* ctx = EVP_PKEY_CTX_new(private_key, NULL);
    if (!ctx) throw std::runtime_error("Failed to create decrypt context");

    if (EVP_PKEY_decrypt_init(ctx) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw std::runtime_error("Failed to init decryption");
    }

    if (EVP_PKEY_CTX_set_rsa_padding(ctx, RSA_PKCS1_OAEP_PADDING) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw std::runtime_error("Failed to set padding");
    }

    size_t outlen;
    if (EVP_PKEY_decrypt(ctx, NULL, &outlen, ciphertext.data(), ciphertext.size()) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw std::runtime_error("Failed to determine plaintext length");
    }

    std::vector<uint8_t> plaintext(outlen);
    if (EVP_PKEY_decrypt(ctx, plaintext.data(), &outlen, ciphertext.data(), ciphertext.size()) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw std::runtime_error("Failed to decrypt");
    }

    EVP_PKEY_CTX_free(ctx);

    uint64_t value = 0;
    memcpy(&value, plaintext.data(), std::min(outlen, sizeof(uint64_t)));
    return value;
}

void generate_test_data(SharedMemory& mem, uint64_t num_devices, uint64_t records_per_device) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> value_dist(1, 1000000);

    mem.private_key = generate_rsa_keypair();
    mem.public_key = mem.private_key;

    uint64_t total_records = num_devices * records_per_device;
    mem.all_data.resize(total_records);

    for (uint64_t device = 0; device < num_devices; ++device) {
        for (uint64_t record = 0; record < records_per_device; ++record) {
            EncryptedData& data = mem.all_data[device * records_per_device + record];
            data.device_id = static_cast<uint32_t>(device);

            uint64_t value = value_dist(gen);
            uint8_t plaintext[8];
            memcpy(plaintext, &value, sizeof(value));

            data.encrypted_value = rsa_encrypt(mem.public_key, plaintext, sizeof(value));
        }
    }

    mem.round_sums.reserve(records_per_device);
}

long long get_time_us() {
#ifdef __gramine__
    struct timespec ts;
    gramine_clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000000LL + ts.tv_nsec / 1000LL;
#else
    auto now = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count();
#endif
}

void run_benchmark(const std::string& platform, SharedMemory& mem, uint64_t num_devices, 
                   uint64_t records_per_device) {

    long long setup_start_us = get_time_us();
    long long setup_time_us = get_time_us() - setup_start_us;

    const uint64_t total_records = num_devices * records_per_device;
    std::vector<size_t> access_order(total_records);
    for (size_t i = 0; i < total_records; ++i) access_order[i] = i;
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(access_order.begin(), access_order.end(), g);

    long long benchmark_start_us = get_time_us();

    for (uint64_t round = 0; round < records_per_device; ++round) {
        uint64_t round_sum = 0;

        for (uint64_t device = 0; device < num_devices; ++device) {
            size_t idx = access_order[round * num_devices + device];
            uint64_t decrypted_value = rsa_decrypt(mem.private_key, mem.all_data[idx].encrypted_value);
            round_sum += decrypted_value;
        }

        mem.round_sums.push_back(round_sum);
    }

    long long benchmark_time_us = get_time_us() - benchmark_start_us;

    long long final_start_us = get_time_us();
    volatile uint64_t total_sum = 0;
    for (uint64_t sum : mem.round_sums) {
        total_sum += sum;
    }
    long long final_time_us = get_time_us() - final_start_us;

    double avg_per_round_us = (records_per_device > 0) ? 
        (double)benchmark_time_us / records_per_device : 0.0;
    double avg_final_us = (records_per_device > 0) ? 
        (double)final_time_us / records_per_device : 0.0;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Total sum: " << total_sum << std::endl;
    std::cout << "\nSetup (one-time):" << std::endl;
    std::cout << "  Key generation: " << setup_time_us << " us" << std::endl;
    std::cout << "\nPer-round averages:" << std::endl;
    std::cout << "  Decrypt+Sum: " << avg_per_round_us << " us" << std::endl;
    std::cout << "  Final agg:   " << avg_final_us << " us" << std::endl;
    std::cout << "  Total:       " << (avg_per_round_us + avg_final_us) << " us" << std::endl;

    const std::string results_filename = "rsa_results.csv";
    std::ofstream results_file;
    std::ifstream check_file(results_filename);
    bool file_exists = check_file.good();
    check_file.close();
    results_file.open(results_filename, std::ios_base::app);

    if (!file_exists) {
        results_file << "platform,num_devices,records_per_device,"
                     << "setup_us,avg_decrypt_sum_per_round_us,avg_final_per_round_us,avg_total_per_round_us\n";
    }

    results_file << platform << "," << num_devices << "," << records_per_device << "," 
                 << setup_time_us << "," << avg_per_round_us << "," 
                 << avg_final_us << "," << (avg_per_round_us + avg_final_us) << "\n";
    results_file.close();
}

int main(int argc, char* argv[]) {
    setenv("MALLOC_ARENA_MAX", "1", 1);

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <command> [args...]" << std::endl;
        return 1;
    }

    std::string command = argv[1];
    if (command != "generate" && command != "native" && command != "sgx") {
        std::cerr << "Unknown command: " << command << std::endl;
        return 1;
    }

    if (argc != 4) {
        std::cerr << "Error: Command requires <num_devices> and <records_per_device>." << std::endl;
        return 1;
    }

    uint64_t num_devices = std::stoull(argv[2]);
    uint64_t records_per_device = std::stoull(argv[3]);

    SharedMemory mem;

    if (command == "generate") {
        std::cout << "Generating test data..." << std::endl;
        generate_test_data(mem, num_devices, records_per_device);
        std::cout << "Data generation complete." << std::endl;
        EVP_PKEY_free(mem.private_key);
        return 0;
    }

    std::cout << "Generating and loading data into memory..." << std::endl;
    generate_test_data(mem, num_devices, records_per_device);
    std::cout << "Memory initialized. Starting benchmark..." << std::endl;

    run_benchmark(command, mem, num_devices, records_per_device);

    EVP_PKEY_free(mem.private_key);

    return 0;
}
