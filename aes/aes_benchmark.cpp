#ifdef __gramine__
#include <gramine/api.hh>
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <random>
#include <cstring>
#include <openssl/rand.h>
#include <openssl/evp.h>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <cstdlib>
#include <algorithm>

struct EncryptedData {
    uint32_t device_id;
    uint8_t encrypted_value[16];
    uint8_t iv[16];
};

struct SharedMemory {
    std::vector<uint8_t> encrypted_keys;
    std::vector<EncryptedData> all_data;
    std::vector<uint64_t> round_sums;
};

std::vector<uint8_t> encrypt_with_gcm(const std::vector<uint8_t>& plaintext, const uint8_t* master_key) {
    EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
    if (!ctx) throw std::runtime_error("Failed to create EVP context");
    uint8_t iv[12];
    RAND_bytes(iv, 12);
    if (EVP_EncryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, master_key, iv) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        throw std::runtime_error("Failed to init encryption");
    }
    std::vector<uint8_t> ciphertext(plaintext.size() + 16);
    int len;
    if (EVP_EncryptUpdate(ctx, ciphertext.data(), &len, plaintext.data(), plaintext.size()) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        throw std::runtime_error("Failed to encrypt");
    }
    int ciphertext_len = len;
    if (EVP_EncryptFinal_ex(ctx, ciphertext.data() + len, &len) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        throw std::runtime_error("Failed to finalize encryption");
    }
    ciphertext_len += len;
    uint8_t tag[16];
    EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_GET_TAG, 16, tag);
    EVP_CIPHER_CTX_free(ctx);
    std::vector<uint8_t> output(12 + ciphertext_len + 16);
    memcpy(output.data(), iv, 12);
    memcpy(output.data() + 12, ciphertext.data(), ciphertext_len);
    memcpy(output.data() + 12 + ciphertext_len, tag, 16);
    return output;
}

std::vector<uint8_t> decrypt_with_gcm(const std::vector<uint8_t>& encrypted_data, const uint8_t* master_key) {
    if (encrypted_data.size() < 12 + 16) throw std::runtime_error("Invalid encrypted data size");
    const uint8_t* iv = encrypted_data.data();
    const uint8_t* ciphertext = encrypted_data.data() + 12;
    const uint8_t* tag = encrypted_data.data() + encrypted_data.size() - 16;
    int ciphertext_len = encrypted_data.size() - 12 - 16;
    EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
    if (!ctx) throw std::runtime_error("Failed to create EVP context");
    if (EVP_DecryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, master_key, iv) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        throw std::runtime_error("Failed to init decryption");
    }
    std::vector<uint8_t> plaintext(ciphertext_len);
    int len;
    if (EVP_DecryptUpdate(ctx, plaintext.data(), &len, ciphertext, ciphertext_len) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        throw std::runtime_error("Failed to decrypt");
    }
    int plaintext_len = len;
    EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_TAG, 16, const_cast<uint8_t*>(tag));
    if (EVP_DecryptFinal_ex(ctx, plaintext.data() + len, &len) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        throw std::runtime_error("Failed to finalize decryption (tag mismatch?)");
    }
    plaintext_len += len;
    EVP_CIPHER_CTX_free(ctx);
    plaintext.resize(plaintext_len);
    return plaintext;
}

class AESDecryptor {
private:
    std::unordered_map<uint32_t, std::vector<uint8_t>> device_keys;
    EVP_CIPHER_CTX* reusable_ctx;
public:
    AESDecryptor() {
        reusable_ctx = EVP_CIPHER_CTX_new();
        if (!reusable_ctx) throw std::runtime_error("Failed to create EVP context");
    }
    ~AESDecryptor() {
        if (reusable_ctx) EVP_CIPHER_CTX_free(reusable_ctx);
    }
    void add_device_key(uint32_t device_id, const std::vector<uint8_t>& key) {
        device_keys[device_id] = key;
    }
    uint64_t decrypt_value(const EncryptedData& data) {
        auto key_it = device_keys.find(data.device_id);
        if (key_it == device_keys.end()) {
            throw std::runtime_error("Unknown device ID");
        }
        EVP_CIPHER_CTX_reset(reusable_ctx);
        if (EVP_DecryptInit_ex(reusable_ctx, EVP_aes_128_cbc(), NULL,
                              key_it->second.data(), data.iv) != 1) {
            throw std::runtime_error("Failed to init decryption");
        }
        EVP_CIPHER_CTX_set_padding(reusable_ctx, 0);
        uint8_t decrypted[16];
        int len;
        if (EVP_DecryptUpdate(reusable_ctx, decrypted, &len, data.encrypted_value, 16) != 1) {
            throw std::runtime_error("Failed to decrypt update");
        }
        if (EVP_DecryptFinal_ex(reusable_ctx, decrypted + len, &len) != 1) {
            throw std::runtime_error("Failed to finalize decryption");
        }
        return *reinterpret_cast<uint64_t*>(decrypted);
    }
};

void generate_test_data(const std::string& data_filename, const std::string& key_filename, uint64_t num_devices, uint64_t records_per_device, const uint8_t* master_key) {
    std::ofstream data_file(data_filename, std::ios::binary);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> value_dist(1, 1000000);
    std::vector<std::vector<uint8_t>> device_keys(num_devices, std::vector<uint8_t>(16));
    for (uint64_t i = 0; i < num_devices; ++i) {
        RAND_bytes(device_keys[i].data(), 16);
    }
    std::vector<uint8_t> key_data(sizeof(uint32_t) + num_devices * 16);
    uint32_t num_keys = static_cast<uint32_t>(num_devices);
    memcpy(key_data.data(), &num_keys, sizeof(uint32_t));
    for (uint64_t i = 0; i < num_devices; ++i) {
        memcpy(key_data.data() + sizeof(uint32_t) + i * 16, device_keys[i].data(), 16);
    }
    auto encrypted_keys = encrypt_with_gcm(key_data, master_key);
    std::ofstream key_file(key_filename, std::ios::binary);
    key_file.write(reinterpret_cast<const char*>(encrypted_keys.data()), encrypted_keys.size());
    uint64_t total_records = num_devices * records_per_device;
    data_file.write(reinterpret_cast<const char*>(&total_records), sizeof(total_records));
    for (uint64_t device = 0; device < num_devices; ++device) {
        for (uint64_t record = 0; record < records_per_device; ++record) {
            EncryptedData data;
            data.device_id = static_cast<uint32_t>(device);
            RAND_bytes(data.iv, 16);
            uint8_t plaintext[16] = {0};
            uint64_t value = value_dist(gen);
            memcpy(plaintext, &value, sizeof(value));
            EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
            if (!ctx) throw std::runtime_error("Failed to create EVP context");
            if (EVP_EncryptInit_ex(ctx, EVP_aes_128_cbc(), NULL, device_keys[device].data(), data.iv) != 1) {
                EVP_CIPHER_CTX_free(ctx);
                throw std::runtime_error("Failed to init encryption");
            }
            EVP_CIPHER_CTX_set_padding(ctx, 0);
            int len;
            if (EVP_EncryptUpdate(ctx, data.encrypted_value, &len, plaintext, 16) != 1) {
                EVP_CIPHER_CTX_free(ctx);
                throw std::runtime_error("Failed to encrypt update");
            }
            if (EVP_EncryptFinal_ex(ctx, data.encrypted_value + len, &len) != 1) {
                EVP_CIPHER_CTX_free(ctx);
                throw std::runtime_error("Failed to finalize encryption");
            }
            EVP_CIPHER_CTX_free(ctx);
            data_file.write(reinterpret_cast<const char*>(&data), sizeof(data));
        }
    }
}

void initialize_shared_memory(SharedMemory& mem, const std::string& data_filename, const std::string& key_filename, uint64_t num_devices, uint64_t records_per_device) {
    std::ifstream key_ifile(key_filename, std::ios::binary);
    if (!key_ifile) throw std::runtime_error("Key file not found. Run 'generate' first.");
    key_ifile.seekg(0, std::ios::end);
    size_t key_file_size = key_ifile.tellg();
    key_ifile.seekg(0, std::ios::beg);
    mem.encrypted_keys.resize(key_file_size);
    key_ifile.read(reinterpret_cast<char*>(mem.encrypted_keys.data()), key_file_size);
    key_ifile.close();

    std::ifstream data_ifile(data_filename, std::ios::binary);
    if (!data_ifile) throw std::runtime_error("Data file not found. Run 'generate' first.");
    uint64_t file_total_records;
    data_ifile.read(reinterpret_cast<char*>(&file_total_records), sizeof(file_total_records));
    const uint64_t total_records = num_devices * records_per_device;
    if (file_total_records != total_records) {
        throw std::runtime_error("Mismatched record count in data file.");
    }
    mem.all_data.resize(total_records);
    data_ifile.read(reinterpret_cast<char*>(mem.all_data.data()), total_records * sizeof(EncryptedData));
    data_ifile.close();

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
                   uint64_t records_per_device, const uint8_t* master_key) {

    // ============ SETUP PHASE (ONE-TIME, REPORTED SEPARATELY) ============
    long long setup_start_us = get_time_us();

    std::vector<uint8_t> decrypted_key_data = decrypt_with_gcm(mem.encrypted_keys, master_key);

    uint32_t num_keys;
    memcpy(&num_keys, decrypted_key_data.data(), sizeof(uint32_t));
    AESDecryptor decryptor;
    for (uint32_t i = 0; i < num_keys; ++i) {
        std::vector<uint8_t> key(16);
        memcpy(key.data(), decrypted_key_data.data() + sizeof(uint32_t) + i * 16, 16);
        decryptor.add_device_key(i, key);
    }

    long long setup_time_us = get_time_us() - setup_start_us;

    // ============ PREPARE ACCESS PATTERN ============
    const uint64_t total_records = num_devices * records_per_device;
    std::vector<size_t> access_order(total_records);
    for (size_t i = 0; i < total_records; ++i) access_order[i] = i;
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(access_order.begin(), access_order.end(), g);

    // ============ BENCHMARK PHASE (PER-ROUND TIMING) ============
    long long benchmark_start_us = get_time_us();

    for (uint64_t round = 0; round < records_per_device; ++round) {
        uint64_t round_sum = 0;

        for (uint64_t device = 0; device < num_devices; ++device) {
            size_t idx = access_order[round * num_devices + device];
            uint64_t decrypted_value = decryptor.decrypt_value(mem.all_data[idx]);
            round_sum += decrypted_value;
        }

        mem.round_sums.push_back(round_sum);
    }

    long long benchmark_time_us = get_time_us() - benchmark_start_us;

    // ============ FINAL AGGREGATION ============
    long long final_start_us = get_time_us();
    volatile uint64_t total_sum = 0;
    for (uint64_t sum : mem.round_sums) {
        total_sum += sum;
    }
    long long final_time_us = get_time_us() - final_start_us;

    // ============ COMPUTE METRICS ============
    double avg_per_round_us = (records_per_device > 0) ? 
        (double)benchmark_time_us / records_per_device : 0.0;
    double avg_final_us = (records_per_device > 0) ? 
        (double)final_time_us / records_per_device : 0.0;

    // ============ REPORT RESULTS ============
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Total sum: " << total_sum << std::endl;
    std::cout << "\nSetup (one-time):" << std::endl;
    std::cout << "  Key decryption: " << setup_time_us << " us" << std::endl;
    std::cout << "\nPer-round averages:" << std::endl;
    std::cout << "  Decrypt+Sum: " << avg_per_round_us << " us" << std::endl;
    std::cout << "  Final agg:   " << avg_final_us << " us" << std::endl;
    std::cout << "  Total:       " << (avg_per_round_us + avg_final_us) << " us" << std::endl;

    // ============ SAVE TO CSV ============
    const std::string results_filename = "aes_results.csv";
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

    const std::string data_filename = "aes_encrypted_data.bin";
    const std::string key_filename = "aes_encrypted_keys.bin";
    uint8_t master_key[32] = {0};

    if (command == "generate") {
        std::cout << "Generating test data..." << std::endl;
        generate_test_data(data_filename, key_filename, num_devices, records_per_device, master_key);
        std::cout << "Data generation complete." << std::endl;
        return 0;
    }

    std::cout << "Initializing shared memory..." << std::endl;
    SharedMemory mem;
    try {
        initialize_shared_memory(mem, data_filename, key_filename, num_devices, records_per_device);
    } catch (const std::exception& e) {
        std::cerr << "Error initializing memory: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Memory initialized. Starting benchmark..." << std::endl;

    run_benchmark(command, mem, num_devices, records_per_device, master_key);

    return 0;
}
