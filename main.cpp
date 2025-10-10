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

// --- UNCHANGED FUNCTIONS (encrypt_with_gcm, decrypt_with_gcm, etc.) ---
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

struct EncryptedData {
    uint32_t device_id;
    uint8_t encrypted_value[16];
    uint8_t iv[16];
};

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

void generate_terse_data(const std::string& data_filename, uint64_t num_devices, uint64_t records_per_device) {
    std::ofstream data_file(data_filename, std::ios::binary);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> value_dist(1, 1000000);
    uint64_t total_records = num_devices * records_per_device;
    data_file.write(reinterpret_cast<const char*>(&total_records), sizeof(total_records));
    for (uint64_t i = 0; i < total_records; ++i) {
        uint64_t value = value_dist(gen);
        data_file.write(reinterpret_cast<const char*>(&value), sizeof(value));
    }
}

int main(int argc, char* argv[]) {
    setenv("MALLOC_ARENA_MAX", "1", 1);

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <command> [args...]" << std::endl;
        std::cerr << "Commands:" << std::endl;
        std::cerr << "  generate <num_devices> <records_per_device>" << std::endl;
        std::cerr << "  generate-terse <num_devices> <records_per_device>" << std::endl;
        std::cerr << "  <platform> <num_devices> <records_per_device> (e.g., native, sgx, terse-native, terse-sgx)" << std::endl;
        std::cerr << "  terse-hybrid-native-sum-looped <num_devices> <records_per_device>" << std::endl;
        std::cerr << "  terse-hybrid-sgx-lookup-looped <num_devices> <records_per_device>" << std::endl;
        return 1;
    }

    std::string command = argv[1];
    uint64_t num_devices = 0;
    uint64_t records_per_device = 0;

    if (command == "generate" || command == "generate-terse" || command == "native" || command == "sgx" || command == "terse-native" || command == "terse-sgx" || command == "terse-hybrid-native-sum-looped" || command == "terse-hybrid-sgx-lookup-looped") {
        if (argc != 4) {
            std::cerr << "Error: Command '" << command << "' requires <num_devices> and <records_per_device>." << std::endl;
            return 1;
        }
        try {
            num_devices = std::stoull(argv[2]);
            records_per_device = std::stoull(argv[3]);
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid number format for devices or records. " << e.what() << std::endl;
            return 1;
        }
    }

    const std::string data_file = "encrypted_data.bin";
    const std::string key_file = "encrypted_keys.bin";
    const std::string terse_data_file = "terse_data.bin";
    const uint64_t total_records = num_devices * records_per_device;

    uint8_t master_key[32] = {0}; // Master key (unchanged)

    if (command == "generate") {
        generate_test_data(data_file, key_file, num_devices, records_per_device, master_key);
        return 0;
    }
    if (command == "generate-terse") {
        generate_terse_data(terse_data_file, num_devices, records_per_device);
        return 0;
    }

    std::string platform = command;
    const size_t LOOKUP_TABLE_SIZE = 16000;
    std::vector<uint64_t> lookup_table(LOOKUP_TABLE_SIZE);
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<uint64_t> value_dist(1, 1000);
        for (size_t i = 0; i < LOOKUP_TABLE_SIZE; ++i) {
            lookup_table[i] = value_dist(gen);
        }
    }

    // --- NEW, ROBUST HYBRID LOGIC ---
    // These modes are special: they time themselves and exit immediately.
    if (platform == "terse-hybrid-native-sum-looped" || platform == "terse-hybrid-sgx-lookup-looped") {
        long long duration_us;
        volatile uint64_t sum = 0; // Use volatile to prevent optimization

        #ifdef __gramine__
            struct timespec start_ts, end_ts;
            gramine_clock_gettime(CLOCK_MONOTONIC, &start_ts);
        #else
            auto start = std::chrono::high_resolution_clock::now();
        #endif

        if (platform == "terse-hybrid-native-sum-looped") {
            std::ifstream data_ifile(terse_data_file, std::ios::binary);
            if (!data_ifile) { std::cerr << "Terse data file not found." << std::endl; return 1; }
            uint64_t file_total_records;
            data_ifile.read(reinterpret_cast<char*>(&file_total_records), sizeof(file_total_records));
            const size_t BATCH_SIZE = 1000;
            std::vector<uint64_t> batch(BATCH_SIZE);
            for (uint64_t round = 0; round < records_per_device; ++round) {
                uint64_t round_sum = 0;
                uint64_t devices_processed = 0;
                while (devices_processed < num_devices) {
                    size_t devices_to_read = std::min(BATCH_SIZE, num_devices - devices_processed);
                    data_ifile.read(reinterpret_cast<char*>(batch.data()), devices_to_read * sizeof(uint64_t));
                    for (size_t i = 0; i < devices_to_read; ++i) {
                        round_sum += batch[i];
                    }
                    devices_processed += devices_to_read;
                }
                sum += round_sum;
            }
        } else { // terse-hybrid-sgx-lookup-looped
            volatile uint64_t dummy_sum = 12345;
            for (uint64_t round = 0; round < records_per_device; ++round) {
                uint64_t lookup_value = lookup_table[dummy_sum % LOOKUP_TABLE_SIZE];
                sum += lookup_value;
            }
        }

        #ifdef __gramine__
            gramine_clock_gettime(CLOCK_MONOTONIC, &end_ts);
            duration_us = (end_ts.tv_sec - start_ts.tv_sec) * 1000000LL + (end_ts.tv_nsec - start_ts.tv_nsec) / 1000LL;
        #else
            auto end = std::chrono::high_resolution_clock::now();
            duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        #endif

        std::cout << "duration_us:" << duration_us << std::endl;
        return 0;
    }

    // --- Standard Benchmark Logic (for all other platforms) ---
    long long duration_us;
    uint64_t sum = 0;
    double average_round_time_us = 0.0;

    #ifdef __gramine__
        struct timespec start_ts, end_ts;
        gramine_clock_gettime(CLOCK_MONOTONIC, &start_ts);
    #else
        auto start = std::chrono::high_resolution_clock::now();
    #endif

    if (platform == "native" || platform == "sgx") {
        std::ifstream data_ifile(data_file, std::ios::binary);
        if (!data_ifile) { std::cerr << "Run with 'generate' first" << std::endl; return 1; }
        std::ifstream key_ifile(key_file, std::ios::binary);
        if (!key_ifile) { std::cerr << "Key file not found" << std::endl; return 1; }
        key_ifile.seekg(0, std::ios::end);
        size_t key_file_size = key_ifile.tellg();
        key_ifile.seekg(0, std::ios::beg);
        std::vector<uint8_t> encrypted_keys(key_file_size);
        key_ifile.read(reinterpret_cast<char*>(encrypted_keys.data()), key_file_size);
        auto decrypted_key_data = decrypt_with_gcm(encrypted_keys, master_key);
        uint32_t num_keys;
        memcpy(&num_keys, decrypted_key_data.data(), sizeof(uint32_t));
        AESDecryptor decryptor;
        for (uint32_t i = 0; i < num_keys; ++i) {
            std::vector<uint8_t> key(16);
            memcpy(key.data(), decrypted_key_data.data() + sizeof(uint32_t) + i * 16, 16);
            decryptor.add_device_key(i, key);
        }
        uint64_t file_total_records;
        data_ifile.read(reinterpret_cast<char*>(&file_total_records), sizeof(file_total_records));
        long long total_round_time = 0;
        long long round_count = 0;
        uint64_t records_processed = 0;
        const size_t BATCH_SIZE = 1000;
        std::vector<EncryptedData> batch(BATCH_SIZE);
        #ifdef __gramine__
            struct timespec round_start_ts, round_end_ts;
            gramine_clock_gettime(CLOCK_MONOTONIC, &round_start_ts);
        #else
            auto round_start = std::chrono::high_resolution_clock::now();
        #endif
        while (records_processed < total_records) {
            size_t records_to_read = std::min(BATCH_SIZE, total_records - records_processed);
            data_ifile.read(reinterpret_cast<char*>(batch.data()), records_to_read * sizeof(EncryptedData));
            for (size_t i = 0; i < records_to_read; ++i) {
                sum += decryptor.decrypt_value(batch[i]);
                records_processed++;
                if (records_processed % num_devices == 0 && num_devices > 0) {
                    #ifdef __gramine__
                        gramine_clock_gettime(CLOCK_MONOTONIC, &round_end_ts);
                        long long round_duration_us = (round_end_ts.tv_sec - round_start_ts.tv_sec) * 1000000LL + (round_end_ts.tv_nsec - round_start_ts.tv_nsec) / 1000LL;
                        round_start_ts = round_end_ts;
                    #else
                        auto round_end = std::chrono::high_resolution_clock::now();
                        long long round_duration_us = std::chrono::duration_cast<std::chrono::microseconds>(round_end - round_start).count();
                        round_start = round_end;
                    #endif
                    total_round_time += round_duration_us;
                    round_count++;
                }
            }
        }
        if (round_count > 0) average_round_time_us = static_cast<double>(total_round_time) / round_count;
    } else if (platform == "terse-native" || platform == "terse-sgx") {
        std::ifstream data_ifile(terse_data_file, std::ios::binary);
        if (!data_ifile) { std::cerr << "Run with 'generate-terse' first" << std::endl; return 1; }
        uint64_t file_total_records;
        data_ifile.read(reinterpret_cast<char*>(&file_total_records), sizeof(file_total_records));
        long long total_round_time = 0;
        long long round_count = 0;
        const size_t BATCH_SIZE = 1000;
        std::vector<uint64_t> batch(BATCH_SIZE);
        for (uint64_t round = 0; round < records_per_device; ++round) {
            #ifdef __gramine__
                struct timespec round_start_ts, round_end_ts;
                gramine_clock_gettime(CLOCK_MONOTONIC, &round_start_ts);
            #else
                auto round_start = std::chrono::high_resolution_clock::now();
            #endif
            uint64_t round_sum = 0;
            uint64_t devices_processed = 0;
            while (devices_processed < num_devices) {
                size_t devices_to_read = std::min(BATCH_SIZE, num_devices - devices_processed);
                data_ifile.read(reinterpret_cast<char*>(batch.data()), devices_to_read * sizeof(uint64_t));
                for (size_t i = 0; i < devices_to_read; ++i) round_sum += batch[i];
                devices_processed += devices_to_read;
            }
            uint64_t lookup_value = lookup_table[round_sum % LOOKUP_TABLE_SIZE];
            sum += round_sum + lookup_value;
            #ifdef __gramine__
                gramine_clock_gettime(CLOCK_MONOTONIC, &round_end_ts);
                long long round_duration_us = (round_end_ts.tv_sec - round_start_ts.tv_sec) * 1000000LL + (round_end_ts.tv_nsec - round_start_ts.tv_nsec) / 1000LL;
            #else
                auto round_end = std::chrono::high_resolution_clock::now();
                long long round_duration_us = std::chrono::duration_cast<std::chrono::microseconds>(round_end - round_start).count();
            #endif
            total_round_time += round_duration_us;
            round_count++;
        }
        if (round_count > 0) average_round_time_us = static_cast<double>(total_round_time) / round_count;
    } else {
        std::cerr << "Unknown platform: " << platform << std::endl;
        return 1;
    }

    #ifdef __gramine__
        gramine_clock_gettime(CLOCK_MONOTONIC, &end_ts);
        duration_us = (end_ts.tv_sec - start_ts.tv_sec) * 1000000LL + (end_ts.tv_nsec - start_ts.tv_nsec) / 1000LL;
    #else
        auto end = std::chrono::high_resolution_clock::now();
        duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    #endif

    std::cout << "Total sum: " << sum << std::endl;
    std::cout << "Time: " << duration_us << " microseconds" << std::endl;
    if (average_round_time_us > 0) {
        std::cout << "Average round time: " << std::fixed << std::setprecision(2) << average_round_time_us << " microseconds" << std::endl;
    }
    double ops_per_sec = (duration_us > 0) ? (total_records * 1000000.0 / duration_us) : 0.0;
    std::cout << "Rate: " << std::fixed << std::setprecision(2) << ops_per_sec << " ops/sec" << std::endl;
    const std::string results_filename = "results.csv";
    std::ofstream results_file;
    std::ifstream check_file(results_filename);
    bool file_exists = check_file.good();
    results_file.open(results_filename, std::ios_base::app);
    if (!file_exists) {
        results_file << "platform,num_devices,records_per_device,duration_us,ops_per_sec,avg_round_time_us\n";
    }
    results_file << platform << "," << num_devices << "," << records_per_device << "," << duration_us << ","
                 << std::fixed << std::setprecision(2) << ops_per_sec << "," << average_round_time_us << "\n";
    results_file.close();
    std::cout << "Results appended to " << results_filename << std::endl;

    return 0;
}
