#ifdef __gramine__
#include <gramine/api.hh>
#endif

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <cstring>
#include <openssl/ec.h>
#include <openssl/evp.h>
#include <openssl/err.h>
#include <openssl/rand.h>
#include <openssl/kdf.h>
#include <openssl/x509.h>
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

EVP_PKEY* generate_ec_keypair() {
    EVP_PKEY* pkey = NULL;
    EVP_PKEY_CTX* pctx = EVP_PKEY_CTX_new_id(EVP_PKEY_EC, NULL);
    if (!pctx) throw std::runtime_error("Failed to create EC context");

    if (EVP_PKEY_keygen_init(pctx) <= 0) {
        EVP_PKEY_CTX_free(pctx);
        throw std::runtime_error("Failed to init keygen");
    }

    if (EVP_PKEY_CTX_set_ec_paramgen_curve_nid(pctx, NID_X9_62_prime256v1) <= 0) {
        EVP_PKEY_CTX_free(pctx);
        throw std::runtime_error("Failed to set curve");
    }

    if (EVP_PKEY_keygen(pctx, &pkey) <= 0) {
        EVP_PKEY_CTX_free(pctx);
        throw std::runtime_error("Failed to generate keypair");
    }

    EVP_PKEY_CTX_free(pctx);
    return pkey;
}

std::vector<uint8_t> derive_shared_secret(EVP_PKEY* private_key, EVP_PKEY* peer_public_key) {
    EVP_PKEY_CTX* ctx = EVP_PKEY_CTX_new(private_key, NULL);
    if (!ctx) throw std::runtime_error("Failed to create derive context");

    if (EVP_PKEY_derive_init(ctx) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw std::runtime_error("Failed to init derive");
    }

    if (EVP_PKEY_derive_set_peer(ctx, peer_public_key) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw std::runtime_error("Failed to set peer");
    }

    size_t secret_len = 0;
    if (EVP_PKEY_derive(ctx, NULL, &secret_len) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw std::runtime_error("Failed to determine secret length");
    }

    std::vector<uint8_t> shared_secret(secret_len);
    if (EVP_PKEY_derive(ctx, shared_secret.data(), &secret_len) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw std::runtime_error("Failed to derive shared secret");
    }
    EVP_PKEY_CTX_free(ctx);

    std::vector<uint8_t> key_material(32);
    unsigned int md_len = 32;
    if (!EVP_Digest(shared_secret.data(), shared_secret.size(), 
                    key_material.data(), &md_len, EVP_sha256(), NULL)) {
        throw std::runtime_error("Failed to hash shared secret");
    }

    return key_material;
}

std::vector<uint8_t> ec_public_key_to_bytes(EVP_PKEY* pkey) {
    size_t len = 0;
    if (EVP_PKEY_get_raw_public_key(pkey, NULL, &len) <= 0) {
        unsigned char* buf = NULL;
        int buf_len = i2d_PUBKEY(pkey, &buf);
        if (buf_len <= 0) throw std::runtime_error("Failed to serialize public key");
        std::vector<uint8_t> result(buf, buf + buf_len);
        OPENSSL_free(buf);
        return result;
    }

    std::vector<uint8_t> result(len);
    if (EVP_PKEY_get_raw_public_key(pkey, result.data(), &len) <= 0) {
        throw std::runtime_error("Failed to get raw public key");
    }
    return result;
}

EVP_PKEY* ec_public_key_from_bytes(const std::vector<uint8_t>& bytes) {
    const unsigned char* p = bytes.data();
    EVP_PKEY* pkey = d2i_PUBKEY(NULL, &p, bytes.size());
    if (!pkey) throw std::runtime_error("Failed to deserialize public key");
    return pkey;
}

std::vector<uint8_t> ecies_encrypt(EVP_PKEY* public_key, const uint8_t* plaintext, size_t plaintext_len) {
    EVP_PKEY* ephemeral_key = generate_ec_keypair();
    std::vector<uint8_t> aes_key = derive_shared_secret(ephemeral_key, public_key);

    uint8_t iv[12];
    if (RAND_bytes(iv, 12) != 1) {
        EVP_PKEY_free(ephemeral_key);
        throw std::runtime_error("Failed to generate IV");
    }

    EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
    if (!ctx) {
        EVP_PKEY_free(ephemeral_key);
        throw std::runtime_error("Failed to create cipher context");
    }

    if (EVP_EncryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, aes_key.data(), iv) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        EVP_PKEY_free(ephemeral_key);
        throw std::runtime_error("Failed to init encryption");
    }

    std::vector<uint8_t> ciphertext(plaintext_len + EVP_CIPHER_block_size(EVP_aes_256_gcm()));
    int len = 0;
    if (EVP_EncryptUpdate(ctx, ciphertext.data(), &len, plaintext, plaintext_len) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        EVP_PKEY_free(ephemeral_key);
        throw std::runtime_error("Failed to encrypt");
    }
    int ciphertext_len = len;

    if (EVP_EncryptFinal_ex(ctx, ciphertext.data() + len, &len) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        EVP_PKEY_free(ephemeral_key);
        throw std::runtime_error("Failed to finalize encryption");
    }
    ciphertext_len += len;
    ciphertext.resize(ciphertext_len);

    uint8_t tag[16];
    if (EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_GET_TAG, 16, tag) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        EVP_PKEY_free(ephemeral_key);
        throw std::runtime_error("Failed to get tag");
    }
    EVP_CIPHER_CTX_free(ctx);

    std::vector<uint8_t> ephemeral_pub = ec_public_key_to_bytes(ephemeral_key);
    EVP_PKEY_free(ephemeral_key);

    std::vector<uint8_t> output;
    uint32_t pub_len = static_cast<uint32_t>(ephemeral_pub.size());
    output.insert(output.end(), reinterpret_cast<uint8_t*>(&pub_len), 
                  reinterpret_cast<uint8_t*>(&pub_len) + sizeof(pub_len));
    output.insert(output.end(), ephemeral_pub.begin(), ephemeral_pub.end());
    output.insert(output.end(), iv, iv + 12);
    output.insert(output.end(), ciphertext.begin(), ciphertext.end());
    output.insert(output.end(), tag, tag + 16);

    return output;
}

uint64_t ecies_decrypt(EVP_PKEY* private_key, const std::vector<uint8_t>& encrypted_data) {
    if (encrypted_data.size() < sizeof(uint32_t) + 12 + 16) {
        throw std::runtime_error("Invalid encrypted data size");
    }

    size_t offset = 0;
    uint32_t pub_len;
    memcpy(&pub_len, encrypted_data.data() + offset, sizeof(pub_len));
    offset += sizeof(pub_len);

    if (encrypted_data.size() < offset + pub_len + 12 + 16) {
        throw std::runtime_error("Invalid encrypted data structure");
    }

    std::vector<uint8_t> ephemeral_pub(encrypted_data.begin() + offset, 
                                        encrypted_data.begin() + offset + pub_len);
    offset += pub_len;

    EVP_PKEY* ephemeral_key = ec_public_key_from_bytes(ephemeral_pub);
    std::vector<uint8_t> aes_key = derive_shared_secret(private_key, ephemeral_key);
    EVP_PKEY_free(ephemeral_key);

    const uint8_t* iv = encrypted_data.data() + offset;
    offset += 12;

    const uint8_t* ciphertext = encrypted_data.data() + offset;
    size_t ciphertext_len = encrypted_data.size() - offset - 16;
    const uint8_t* tag = encrypted_data.data() + encrypted_data.size() - 16;

    EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
    if (!ctx) throw std::runtime_error("Failed to create cipher context");

    if (EVP_DecryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, aes_key.data(), iv) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        throw std::runtime_error("Failed to init decryption");
    }

    std::vector<uint8_t> plaintext(ciphertext_len + EVP_CIPHER_block_size(EVP_aes_256_gcm()));
    int len = 0;
    if (EVP_DecryptUpdate(ctx, plaintext.data(), &len, ciphertext, ciphertext_len) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        throw std::runtime_error("Failed to decrypt");
    }
    int plaintext_len = len;

    if (EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_TAG, 16, const_cast<uint8_t*>(tag)) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        throw std::runtime_error("Failed to set tag");
    }

    if (EVP_DecryptFinal_ex(ctx, plaintext.data() + len, &len) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        throw std::runtime_error("Failed to finalize decryption");
    }
    plaintext_len += len;
    EVP_CIPHER_CTX_free(ctx);

    uint64_t value = 0;
    memcpy(&value, plaintext.data(), std::min(sizeof(uint64_t), static_cast<size_t>(plaintext_len)));
    return value;
}

void generate_test_data(SharedMemory& mem, uint64_t num_devices, uint64_t records_per_device) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> value_dist(1, 1000000);

    mem.private_key = generate_ec_keypair();
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

            data.encrypted_value = ecies_encrypt(mem.public_key, plaintext, sizeof(value));
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

void run_benchmark(const std::string& platform, SharedMemory& mem, uint64_t num_devices, uint64_t records_per_device) {
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
            uint64_t decrypted_value = ecies_decrypt(mem.private_key, mem.all_data[idx].encrypted_value);
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

    double avg_per_round_us = (records_per_device > 0) ? (double)benchmark_time_us / records_per_device : 0.0;
    double avg_final_us = (records_per_device > 0) ? (double)final_time_us / records_per_device : 0.0;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Total sum: " << total_sum << std::endl;
    std::cout << "\nSetup (one-time):" << std::endl;
    std::cout << "  Key generation: " << setup_time_us << " us" << std::endl;
    std::cout << "\nPer-round averages:" << std::endl;
    std::cout << "  Decrypt+Sum: " << avg_per_round_us << " us" << std::endl;
    std::cout << "  Final agg:   " << avg_final_us << " us" << std::endl;
    std::cout << "  Total:       " << (avg_per_round_us + avg_final_us) << " us" << std::endl;

    const std::string results_filename = "ecc_results.csv";
    std::ofstream results_file;
    std::ifstream check_file(results_filename);
    bool file_exists = check_file.good();
    check_file.close();
    results_file.open(results_filename, std::ios_base::app);

    if (!file_exists) {
        results_file << "platform,num_devices,records_per_device,setup_us,avg_decrypt_sum_per_round_us,avg_final_per_round_us,avg_total_per_round_us\n";
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
