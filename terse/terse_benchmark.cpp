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

struct SharedMemory {
    std::vector<uint64_t> input_data;
    std::vector<uint64_t> intermediate_results;
    std::vector<uint64_t> final_results;
    std::vector<uint64_t> lookup_table;
};

void generate_terse_data(SharedMemory& mem, uint64_t num_devices, uint64_t records_per_device) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> value_dist(1, 1000000);
    uint64_t total_records = num_devices * records_per_device;

    mem.input_data.resize(total_records);
    for (uint64_t i = 0; i < total_records; ++i) {
        mem.input_data[i] = value_dist(gen);
    }
}

void initialize_shared_memory(SharedMemory& mem, uint64_t num_devices, uint64_t records_per_device) {
    const size_t LOOKUP_TABLE_SIZE = 16000;

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
    long long total_start_us = get_time_us();

    long long sum_start_us = get_time_us();
    for (uint64_t round = 0; round < records_per_device; ++round) {
        uint64_t round_sum = 0;
        for (uint64_t device = 0; device < num_devices; ++device) {
            round_sum += mem.input_data[device * records_per_device + round];
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

    double avg_sum_us = (records_per_device > 0) ? (double)sum_time_us / records_per_device : 0.0;
    double avg_lookup_us = (records_per_device > 0) ? (double)lookup_time_us / records_per_device : 0.0;
    double avg_total_us = (records_per_device > 0) ? (double)total_time_us / records_per_device : 0.0;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Total: " << total << std::endl;
    std::cout << "Average time per round:" << std::endl;
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
        results_file << "platform,num_devices,records_per_device,avg_sum_per_round_us,avg_lookup_per_round_us,avg_total_per_round_us\n";
    }
    results_file << platform << "," << num_devices << "," << records_per_device << ","
                 << avg_sum_us << "," << avg_lookup_us << "," << avg_total_us << "\n";
    results_file.close();
}

void run_hybrid_native(SharedMemory& mem, uint64_t num_devices, uint64_t records_per_device) {
    long long sum_start_us = get_time_us();
    for (uint64_t round = 0; round < records_per_device; ++round) {
        uint64_t round_sum = 0;
        for (uint64_t device = 0; device < num_devices; ++device) {
            round_sum += mem.input_data[device * records_per_device + round];
        }
        mem.intermediate_results[round] = round_sum;
    }
    long long sum_time_us = get_time_us() - sum_start_us;

    double avg_sum_us = (records_per_device > 0) ? (double)sum_time_us / records_per_device : 0.0;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Native phase - Average per round:" << std::endl;
    std::cout << "  Sum:  " << avg_sum_us << " us" << std::endl;
}

void run_hybrid_sgx(SharedMemory& mem, uint64_t num_devices, uint64_t records_per_device) {
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

    double avg_lookup_us = (records_per_device > 0) ? (double)lookup_time_us / records_per_device : 0.0;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Total: " << total << std::endl;
    std::cout << "SGX phase - Average per round:" << std::endl;
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
        SharedMemory mem;
        generate_terse_data(mem, num_devices, records_per_device);
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

    if (command == "native" || command == "sgx_only") {
        std::cout << "Generating and loading data into memory..." << std::endl;
        generate_terse_data(mem, num_devices, records_per_device);
        initialize_shared_memory(mem, num_devices, records_per_device);
        std::cout << "Memory initialized. Starting benchmark..." << std::endl;
        run_benchmark(command, mem, num_devices, records_per_device);
        std::cout << "Results appended to terse_results.csv" << std::endl;
    } else if (command == "hybrid_native") {
        std::cout << "Generating and loading data into memory..." << std::endl;
        generate_terse_data(mem, num_devices, records_per_device);
        initialize_shared_memory(mem, num_devices, records_per_device);
        std::cout << "Memory initialized. Starting native aggregation..." << std::endl;
        run_hybrid_native(mem, num_devices, records_per_device);
    } else if (command == "hybrid_sgx") {
        std::cout << "Initializing memory for hybrid (SGX part)..." << std::endl;
        initialize_shared_memory(mem, num_devices, records_per_device);
        std::cout << "Memory initialized. Starting SGX lookup..." << std::endl;
        run_hybrid_sgx(mem, num_devices, records_per_device);
    } else {
        std::cerr << "Unknown command: " << command << std::endl;
        return 1;
    }

    return 0;
}
