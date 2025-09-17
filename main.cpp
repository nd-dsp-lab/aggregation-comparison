#ifdef __gramine__
#include <gramine/api.hh>
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <random>
#include <cstring>  // For memcpy
#include <openssl/rand.h>
#include <openssl/evp.h>  // Using EVP for all crypto operations (OpenSSL 3.0 compatible)
#include <cstdint>  // For uint64_t, uint32_t, uint8_t
#include <string>   // For std::stoull
#include <stdexcept> // For std::invalid_argument

// Function to encrypt data with AES-256-GCM (unchanged)
std::vector<uint8_t> encrypt_with_gcm(const std::vector<uint8_t>& plaintext, const uint8_t* master_key) {
    EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
    if (!ctx) throw std::runtime_error("Failed to create EVP context");

    uint8_t iv[12];
    RAND_bytes(iv, 12);

    if (EVP_EncryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, master_key, iv) != 1) {
        EVP_CIPHER_CTX_free(ctx);
        throw std::runtime_error("Failed to init encryption");
    }

    std::vector<uint8_t> ciphertext(plaintext.size() + 16);  // Room for tag
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

    // Prepare output: IV + ciphertext + tag
    std::vector<uint8_t> output(12 + ciphertext_len + 16);
    memcpy(output.data(), iv, 12);
    memcpy(output.data() + 12, ciphertext.data(), ciphertext_len);
    memcpy(output.data() + 12 + ciphertext_len, tag, 16);

    return output;
}

// Function to decrypt data with AES-256-GCM (unchanged)
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
    uint8_t encrypted_value[16]; // AES block size
    uint8_t iv[16];
};

class AESDecryptor {
private:
    std::unordered_map<uint32_t, std::vector<uint8_t>> device_keys;

public:
    void add_device_key(uint32_t device_id, const std::vector<uint8_t>& key) {
        device_keys[device_id] = key;
    }

    uint64_t decrypt_value(const EncryptedData& data) {
        auto key_it = device_keys.find(data.device_id);
        if (key_it == device_keys.end()) {
            throw std::runtime_error("Unknown device ID");
        }

        EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
        if (!ctx) throw std::runtime_error("Failed to create EVP context");

        if (EVP_DecryptInit_ex(ctx, EVP_aes_128_cbc(), NULL, key_it->second.data(), data.iv) != 1) {
            EVP_CIPHER_CTX_free(ctx);
            throw std::runtime_error("Failed to init decryption");
        }
        EVP_CIPHER_CTX_set_padding(ctx, 0);  // Disable padding

        uint8_t decrypted[16 + 16];  // Extra space (though not needed with no padding)
        int len;
        if (EVP_DecryptUpdate(ctx, decrypted, &len, data.encrypted_value, 16) != 1) {
            EVP_CIPHER_CTX_free(ctx);
            throw std::runtime_error("Failed to decrypt update");
        }
        int decrypted_len = len;

        if (EVP_DecryptFinal_ex(ctx, decrypted + len, &len) != 1) {
            EVP_CIPHER_CTX_free(ctx);
            throw std::runtime_error("Failed to finalize decryption");
        }
        decrypted_len += len;

        EVP_CIPHER_CTX_free(ctx);

        // Extract uint64_t from first 8 bytes (assuming decrypted_len >= 8)
        return *reinterpret_cast<uint64_t*>(decrypted);
    }
};

void generate_test_data(const std::string& data_filename, const std::string& key_filename, uint64_t num_devices, uint64_t records_per_device, const uint8_t* master_key) {
    std::ofstream data_file(data_filename, std::ios::binary);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> value_dist(1, 1000000);

    // Generate keys for each device
    std::vector<std::vector<uint8_t>> device_keys(num_devices, std::vector<uint8_t>(16));
    for (uint64_t i = 0; i < num_devices; ++i) {
        RAND_bytes(device_keys[i].data(), 16);
    }

    // Prepare key data as flat buffer: num_keys (4 bytes) + keys (16 bytes each)
    std::vector<uint8_t> key_data(sizeof(uint32_t) + num_devices * 16);
    uint32_t num_keys = static_cast<uint32_t>(num_devices);  // Assume num_devices fits in uint32_t
    memcpy(key_data.data(), &num_keys, sizeof(uint32_t));
    for (uint64_t i = 0; i < num_devices; ++i) {
        memcpy(key_data.data() + sizeof(uint32_t) + i * 16, device_keys[i].data(), 16);
    }

    // Encrypt key data with master key and write to separate key file
    auto encrypted_keys = encrypt_with_gcm(key_data, master_key);
    std::ofstream key_file(key_filename, std::ios::binary);
    key_file.write(reinterpret_cast<const char*>(encrypted_keys.data()), encrypted_keys.size());

    // Generate and encrypt data (now only in data file)
    uint64_t total_records = num_devices * records_per_device;
    data_file.write(reinterpret_cast<const char*>(&total_records), sizeof(total_records));

    for (uint64_t device = 0; device < num_devices; ++device) {
        for (uint64_t record = 0; record < records_per_device; ++record) {
            EncryptedData data;
            data.device_id = static_cast<uint32_t>(device);

            // Generate random IV
            RAND_bytes(data.iv, 16);

            // Create plaintext with random value
            uint8_t plaintext[16] = {0};
            uint64_t value = value_dist(gen);
            memcpy(plaintext, &value, sizeof(value));

            // Encrypt using EVP (OpenSSL 3.0 compatible)
            EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
            if (!ctx) throw std::runtime_error("Failed to create EVP context");

            if (EVP_EncryptInit_ex(ctx, EVP_aes_128_cbc(), NULL, device_keys[device].data(), data.iv) != 1) {
                EVP_CIPHER_CTX_free(ctx);
                throw std::runtime_error("Failed to init encryption");
            }
            EVP_CIPHER_CTX_set_padding(ctx, 0);  // Disable padding

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

// --- NEW: Function to generate data for TERSE benchmark ---
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
    // --- Argument Parsing ---
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <command> [args...]" << std::endl;
        std::cerr << "Commands:" << std::endl;
        std::cerr << "  generate <num_devices> <records_per_device>" << std::endl;
        std::cerr << "  generate-terse <num_devices> <records_per_device>" << std::endl;
        std::cerr << "  <platform> <num_devices> <records_per_device> (e.g., native, sgx, terse-sgx)" << std::endl;
        return 1;
    }

    std::string command = argv[1];
    uint64_t num_devices = 0;
    uint64_t records_per_device = 0;

    // Commands that require device/record counts
    if (command == "generate" || command == "generate-terse" || command == "native" || command == "sgx" || command == "terse-sgx") {
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

    uint8_t master_key[32] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
                              0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F,
                              0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
                              0x18, 0x19, 0x1A, 0x1B, 0x1C, 0x1D, 0x1E, 0x1F};

    // --- Generation Logic ---
    if (command == "generate") {
        std::cout << "Generating standard encrypted test data..." << std::endl;
        generate_test_data(data_file, key_file, num_devices, records_per_device, master_key);
        std::cout << "Generated " << total_records << " encrypted records ("
                  << num_devices << " devices x " << records_per_device << " records/device)" << std::endl;
        std::cout << "Keys stored securely in " << key_file << std::endl;
        return 0;
    }
    if (command == "generate-terse") {
        std::cout << "Generating TERSE test data..." << std::endl;
        generate_terse_data(terse_data_file, num_devices, records_per_device);
        std::cout << "Generated " << total_records << " integer records for TERSE in " << terse_data_file << std::endl;
        return 0;
    }

    // --- Benchmark Execution Logic ---
    std::string platform = command;
    std::cout << "Processing " << total_records << " values ("
              << num_devices << " devices x " << records_per_device << " records/device) on platform: " << platform << std::endl;

    long long duration_us;
    uint64_t sum = 0;

#ifdef __gramine__
    struct timespec start_ts, end_ts;
    gramine_clock_gettime(CLOCK_MONOTONIC, &start_ts);
#else
    auto start = std::chrono::high_resolution_clock::now();
#endif

    if (platform == "native" || platform == "sgx") {
        // --- Standard Decryption Benchmark ---
        std::ifstream data_ifile(data_file, std::ios::binary);
        if (!data_ifile) { std::cerr << "Run with 'generate " << platform << "' argument first to create test data" << std::endl; return 1; }
        std::ifstream key_ifile(key_file, std::ios::binary);
        if (!key_ifile) { std::cerr << "Key file not found: " << key_file << std::endl; return 1; }

        key_ifile.seekg(0, std::ios::end);
        size_t key_file_size = key_ifile.tellg();
        key_ifile.seekg(0, std::ios::beg);
        std::vector<uint8_t> encrypted_keys(key_file_size);
        key_ifile.read(reinterpret_cast<char*>(encrypted_keys.data()), key_file_size);

        auto decrypted_key_data = decrypt_with_gcm(encrypted_keys, master_key);
        if (decrypted_key_data.size() < sizeof(uint32_t)) { std::cerr << "Invalid decrypted key data" << std::endl; return 1; }
        uint32_t num_keys;
        memcpy(&num_keys, decrypted_key_data.data(), sizeof(uint32_t));
        if (decrypted_key_data.size() != sizeof(uint32_t) + num_keys * 16) { std::cerr << "Decrypted key data size mismatch" << std::endl; return 1; }

        AESDecryptor decryptor;
        for (uint32_t i = 0; i < num_keys; ++i) {
            std::vector<uint8_t> key(16);
            memcpy(key.data(), decrypted_key_data.data() + sizeof(uint32_t) + i * 16, 16);
            decryptor.add_device_key(i, key);
        }

        uint64_t file_total_records;
        data_ifile.read(reinterpret_cast<char*>(&file_total_records), sizeof(file_total_records));

        EncryptedData data;
        for (uint64_t i = 0; i < total_records; ++i) {
            data_ifile.read(reinterpret_cast<char*>(&data), sizeof(data));
            sum += decryptor.decrypt_value(data);
        }

    } else if (platform == "terse-sgx") {
        // --- TERSE Benchmark ---
        std::ifstream data_ifile(terse_data_file, std::ios::binary);
        if (!data_ifile) { std::cerr << "Run with 'generate " << platform << "' argument first to create test data" << std::endl; return 1; }

        uint64_t file_total_records;
        data_ifile.read(reinterpret_cast<char*>(&file_total_records), sizeof(file_total_records));

        // Simulate native aggregation by reading and summing values
        uint64_t value;
        for (uint64_t i = 0; i < total_records; ++i) {
            data_ifile.read(reinterpret_cast<char*>(&value), sizeof(value));
            sum += value;
        }

        // Simulate final in-enclave addition of a preset number
        const uint64_t preset_number = 42;
        sum += preset_number;

        #ifndef __gramine__
        std::cerr << "Warning: Running 'terse-sgx' benchmark natively. The entire operation is untrusted." << std::endl;
        #endif

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
    double ops_per_sec = 0.0;
    if (duration_us > 0) {
        ops_per_sec = (total_records * 1000000.0 / duration_us);
        std::cout << "Rate: " << ops_per_sec << " ops/sec" << std::endl;
    }

    const std::string results_filename = "results.csv";
    std::ofstream results_file;
    std::ifstream check_file(results_filename);
    bool file_exists = check_file.good();
    check_file.close();
    results_file.open(results_filename, std::ios_base::app);
    if (results_file.is_open()) {
        if (!file_exists) {
            results_file << "platform,num_devices,records_per_device,duration_us,ops_per_sec\n";
        }
        results_file << platform << "," << num_devices << "," << records_per_device << "," << duration_us << "," << ops_per_sec << "\n";
        results_file.close();
        std::cout << "Results appended to " << results_filename << std::endl;
    } else {
        std::cerr << "Failed to open " << results_filename << " for writing." << std::endl;
    }

    return 0;
}
