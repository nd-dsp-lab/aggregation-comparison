#ifdef __gramine__
#include <gramine/api.hh>
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>
#include <cstring>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <cstdlib>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

const std::string SHM_NAME = "/terse_intermediate_shm";
const std::string SHM_INPUT_NAME = "/terse_input_shm";

struct SharedMemory {
    std::vector<uint64_t> input_data;
    std::vector<uint64_t> intermediate_results;
    std::vector<uint64_t> final_results;
    std::vector<uint64_t> lookup_table;
    uint64_t* input_shm_ptr;
    size_t input_shm_size;
    int input_shm_fd;
};

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

void initialize_shared_memory(SharedMemory& mem, uint64_t num_devices, uint64_t records_per_device, const std::string& data_file, bool load_input = true, bool use_shm_input = false) {
    const uint64_t total_records = num_devices * records_per_device;
    const size_t LOOKUP_TABLE_SIZE = 16000;

    mem.input_shm_ptr = nullptr;
    mem.input_shm_size = 0;
    mem.input_shm_fd = -1;

    if (load_input) {
        if (use_shm_input) {
            mem.input_shm_size = total_records * sizeof(uint64_t);
            shm_unlink(SHM_INPUT_NAME.c_str());
            mem.input_shm_fd = shm_open(SHM_INPUT_NAME.c_str(), O_CREAT | O_RDWR, 0666);
            if (mem.input_shm_fd == -1) {
                throw std::runtime_error("Failed to create input shared memory");
            }

            if (ftruncate(mem.input_shm_fd, mem.input_shm_size) == -1) {
                close(mem.input_shm_fd);
                throw std::runtime_error("Failed to set input shared memory size");
            }

            mem.input_shm_ptr = static_cast<uint64_t*>(mmap(nullptr, mem.input_shm_size, PROT_READ | PROT_WRITE, MAP_SHARED, mem.input_shm_fd, 0));
            if (mem.input_shm_ptr == MAP_FAILED) {
                close(mem.input_shm_fd);
                throw std::runtime_error("Failed to map input shared memory");
            }

            std::ifstream data_ifile(data_file, std::ios::binary);
            if (!data_ifile) {
                munmap(mem.input_shm_ptr, mem.input_shm_size);
                close(mem.input_shm_fd);
                throw std::runtime_error("Data file not found. Run 'generate' first.");
            }
            uint64_t file_total_records;
            data_ifile.read(reinterpret_cast<char*>(&file_total_records), sizeof(file_total_records));
            data_ifile.read(reinterpret_cast<char*>(mem.input_shm_ptr), total_records * sizeof(uint64_t));
            data_ifile.close();
        } else {
            mem.input_data.resize(total_records);
            std::ifstream data_ifile(data_file, std::ios::binary);
            if (!data_ifile) {
                throw std::runtime_error("Data file not found. Run 'generate' first.");
            }
            uint64_t file_total_records;
            data_ifile.read(reinterpret_cast<char*>(&file_total_records), sizeof(file_total_records));
            data_ifile.read(reinterpret_cast<char*>(mem.input_data.data()), total_records * sizeof(uint64_t));
            data_ifile.close();
        }
    }

    mem.intermediate_results.resize(records_per_device);
    mem.final_results.resize(records_per_device);
    mem.lookup_table.resize(LOOKUP_TABLE_SIZE);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> value_dist(1, 1000);
    for (size_t i = 0; i < LOOKUP_TABLE_SIZE; ++i) {
        mem.lookup_table[i] = value_dist(gen);
    }
}

void cleanup_shared_memory(SharedMemory& mem) {
    if (mem.input_shm_ptr != nullptr) {
        munmap(mem.input_shm_ptr, mem.input_shm_size);
        mem.input_shm_ptr = nullptr;
    }
    if (mem.input_shm_fd != -1) {
        close(mem.input_shm_fd);
        shm_unlink(SHM_INPUT_NAME.c_str());
        mem.input_shm_fd = -1;
    }
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

void run_benchmark(const std::string& platform, SharedMemory& mem, uint64_t num_devices, uint64_t records_per_device, bool use_shm_input = false) {
    long long total_start_us = get_time_us();

    uint64_t* input_ptr = use_shm_input ? mem.input_shm_ptr : mem.input_data.data();

    long long read_start_us = get_time_us();
    volatile uint64_t read_check = input_ptr[0];
    (void)read_check;
    long long read_time_us = get_time_us() - read_start_us;

    long long sum_start_us = get_time_us();
    for (uint64_t round = 0; round < records_per_device; ++round) {
        uint64_t round_sum = 0;
        for (uint64_t device = 0; device < num_devices; ++device) {
            round_sum += input_ptr[device * records_per_device + round];
        }
        mem.intermediate_results[round] = round_sum;
    }
    long long sum_time_us = get_time_us() - sum_start_us;

    long long lookup_start_us = get_time_us();
    volatile uint64_t total = 0;
    for (uint64_t round = 0; round < records_per_device; ++round) {
        uint64_t round_sum = mem.intermediate_results[round];
        uint64_t lookup_value = mem.lookup_table[round_sum % mem.lookup_table.size()];
        uint64_t result = round_sum + lookup_value;
        mem.final_results[round] = result;
        total += result;
    }
    long long lookup_time_us = get_time_us() - lookup_start_us;

    long long total_time_us = get_time_us() - total_start_us;

    double avg_read_us = (records_per_device > 0) ? (double)read_time_us / records_per_device : 0.0;
    double avg_sum_us = (records_per_device > 0) ? (double)sum_time_us / records_per_device : 0.0;
    double avg_lookup_us = (records_per_device > 0) ? (double)lookup_time_us / records_per_device : 0.0;
    double avg_total_us = (records_per_device > 0) ? (double)total_time_us / records_per_device : 0.0;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Total: " << total << std::endl;
    std::cout << "Average time per round:" << std::endl;
    std::cout << "  Read:    " << avg_read_us << " us" << std::endl;
    std::cout << "  Sum:     " << avg_sum_us << " us" << std::endl;
    std::cout << "  Lookup:  " << avg_lookup_us << " us" << std::endl;
    std::cout << "  Total:   " << avg_total_us << " us" << std::endl;

    const std::string results_filename = "terse_results.csv";
    std::ofstream results_file;
    std::ifstream check_file(results_filename);
    bool file_exists = check_file.good();
    check_file.close();
    results_file.open(results_filename, std::ios_base::app);
    if (!file_exists) {
        results_file << "platform,num_devices,records_per_device,avg_read_per_round_us,avg_sum_per_round_us,avg_lookup_per_round_us,avg_total_per_round_us\n";
    }
    results_file << platform << "," << num_devices << "," << records_per_device << ","
                 << avg_read_us << "," << avg_sum_us << "," << avg_lookup_us << "," << avg_total_us << "\n";
    results_file.close();
}

void benchmark_hybrid_native(SharedMemory& mem, uint64_t num_devices, uint64_t records_per_device) {
    size_t shm_size = records_per_device * sizeof(uint64_t);

    shm_unlink(SHM_NAME.c_str());
    int shm_fd = shm_open(SHM_NAME.c_str(), O_CREAT | O_RDWR, 0666);
    if (shm_fd == -1) {
        throw std::runtime_error("Failed to create shared memory");
    }

    if (ftruncate(shm_fd, shm_size) == -1) {
        close(shm_fd);
        throw std::runtime_error("Failed to set shared memory size");
    }

    uint64_t* shm_ptr = static_cast<uint64_t*>(mmap(nullptr, shm_size, PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0));
    if (shm_ptr == MAP_FAILED) {
        close(shm_fd);
        throw std::runtime_error("Failed to map shared memory");
    }

    long long read_start_us = get_time_us();
    volatile uint64_t read_check = mem.input_data[0];
    (void)read_check;
    long long read_time_us = get_time_us() - read_start_us;

    long long sum_start_us = get_time_us();
    for (uint64_t round = 0; round < records_per_device; ++round) {
        uint64_t round_sum = 0;
        for (uint64_t device = 0; device < num_devices; ++device) {
            round_sum += mem.input_data[device * records_per_device + round];
        }
        shm_ptr[round] = round_sum;
    }
    long long sum_time_us = get_time_us() - sum_start_us;

    munmap(shm_ptr, shm_size);
    close(shm_fd);

    double avg_read_us = (records_per_device > 0) ? (double)read_time_us / records_per_device : 0.0;
    double avg_sum_us = (records_per_device > 0) ? (double)sum_time_us / records_per_device : 0.0;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Native phase - Average per round:" << std::endl;
    std::cout << "  Read: " << avg_read_us << " us" << std::endl;
    std::cout << "  Sum:  " << avg_sum_us << " us" << std::endl;
}

void benchmark_hybrid_sgx(SharedMemory& mem, uint64_t num_devices, uint64_t records_per_device) {
    size_t shm_size = records_per_device * sizeof(uint64_t);

    long long read_start_us = get_time_us();
    int shm_fd = shm_open(SHM_NAME.c_str(), O_RDONLY, 0666);
    if (shm_fd == -1) {
        throw std::runtime_error("Failed to open shared memory");
    }

    uint64_t* shm_ptr = static_cast<uint64_t*>(mmap(nullptr, shm_size, PROT_READ, MAP_SHARED, shm_fd, 0));
    if (shm_ptr == MAP_FAILED) {
        close(shm_fd);
        throw std::runtime_error("Failed to map shared memory");
    }
    long long read_time_us = get_time_us() - read_start_us;

    long long lookup_start_us = get_time_us();
    volatile uint64_t total = 0;
    for (uint64_t round = 0; round < records_per_device; ++round) {
        uint64_t round_sum = shm_ptr[round];
        uint64_t lookup_value = mem.lookup_table[round_sum % mem.lookup_table.size()];
        uint64_t result = round_sum + lookup_value;
        mem.final_results[round] = result;
        total += result;
    }
    long long lookup_time_us = get_time_us() - lookup_start_us;

    munmap(shm_ptr, shm_size);
    close(shm_fd);
    shm_unlink(SHM_NAME.c_str());

    double avg_read_us = (records_per_device > 0) ? (double)read_time_us / records_per_device : 0.0;
    double avg_lookup_us = (records_per_device > 0) ? (double)lookup_time_us / records_per_device : 0.0;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Total: " << total << std::endl;
    std::cout << "SGX phase - Average per round:" << std::endl;
    std::cout << "  Read:   " << avg_read_us << " us" << std::endl;
    std::cout << "  Lookup: " << avg_lookup_us << " us" << std::endl;
}

int main(int argc, char* argv[]) {
    setenv("MALLOC_ARENA_MAX", "1", 1);

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <command> [args...]" << std::endl;
        return 1;
    }

    std::string command = argv[1];

    if (command == "generate") {
        if (argc != 4) {
            std::cerr << "Error: generate requires <num_devices> and <records_per_device>." << std::endl;
            return 1;
        }
        uint64_t num_devices = std::stoull(argv[2]);
        uint64_t records_per_device = std::stoull(argv[3]);
        generate_terse_data("terse_data.bin", num_devices, records_per_device);
        std::cout << "Generated data for " << num_devices << " devices, " 
                  << records_per_device << " records each." << std::endl;
        return 0;
    }

    if (argc != 4) {
        std::cerr << "Error: Command requires <num_devices> and <records_per_device>." << std::endl;
        return 1;
    }

    uint64_t num_devices = std::stoull(argv[2]);
    uint64_t records_per_device = std::stoull(argv[3]);

    SharedMemory mem;

    if (command == "native") {
        std::cout << "Initializing memory for native..." << std::endl;
        initialize_shared_memory(mem, num_devices, records_per_device, "terse_data.bin", true, false);
        std::cout << "Memory initialized. Starting benchmark..." << std::endl;
        run_benchmark(command, mem, num_devices, records_per_device, false);
        std::cout << "Results appended to terse_results.csv" << std::endl;
    } else if (command == "sgx") {
        std::cout << "Initializing memory for SGX (using untrusted shared memory for input)..." << std::endl;
        initialize_shared_memory(mem, num_devices, records_per_device, "terse_data.bin", true, true);
        std::cout << "Memory initialized. Starting benchmark..." << std::endl;
        run_benchmark(command, mem, num_devices, records_per_device, true);
        cleanup_shared_memory(mem);
        std::cout << "Results appended to terse_results.csv" << std::endl;
    } else if (command == "hybrid_native") {
        std::cout << "Initializing memory for hybrid (native part)..." << std::endl;
        initialize_shared_memory(mem, num_devices, records_per_device, "terse_data.bin");
        std::cout << "Memory initialized. Starting native aggregation..." << std::endl;
        benchmark_hybrid_native(mem, num_devices, records_per_device);
    } else if (command == "hybrid_sgx") {
        std::cout << "Initializing memory for hybrid (SGX part)..." << std::endl;
        initialize_shared_memory(mem, num_devices, records_per_device, "", false);
        std::cout << "Memory initialized. Starting SGX lookup/write..." << std::endl;
        benchmark_hybrid_sgx(mem, num_devices, records_per_device);
    } else {
        std::cerr << "Unknown command: " << command << std::endl;
        return 1;
    }

    return 0;
}
