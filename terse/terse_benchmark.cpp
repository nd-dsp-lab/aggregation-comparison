#include "terse_benchmark.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <random>

uint64_t TERSEBenchmark::get_time_us() {
    return std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now().time_since_epoch()
    ).count();
}

void TERSEBenchmark::generate_all_test_data() {
    DiscreteLaplacian dl;
    uint64_t max_val = (1ULL << config.plain_bits) - 1;

    all_user_inputs.resize(config.records_per_device);
    for (uint64_t round = 0; round < config.records_per_device; round++) {
        all_user_inputs[round].resize(config.num_devices);
        for (uint64_t device = 0; device < config.num_devices; device++) {
            all_user_inputs[round][device] = dl.uniform_64(max_val);
        }
    }
}

TERSEBenchmark::TERSEBenchmark(const BenchmarkConfig& config) : config(config) {
    std::cout << "Initializing TERSE with:" << std::endl;
    std::cout << "  num_devices: " << config.num_devices << std::endl;
    std::cout << "  plain_bits: " << config.plain_bits << std::endl;
    std::cout << "  scale: " << config.scale << std::endl;

    double beta = 1.0 / config.num_devices;
    std::cout << "  beta: " << beta << std::endl;
    std::cout << "  scheme: " << (config.scheme == NS ? "NS" : "MS") << std::endl;

    std::cout << "Creating Aggregator_RNS..." << std::endl;

    agg = new Aggregator_RNS(
        config.plain_bits,
        config.scale,
        config.num_devices,
        config.scheme,
        beta,
        config.plain_bits
    );

    std::cout << "Aggregator_RNS created successfully" << std::endl;

    auto params_pair = agg->parms_ptrs();
    ctext_parms = params_pair.first;
    plain_parms = params_pair.second;

    agg_key = new Polynomial(ctext_parms);
    agg->secret_keys(*agg_key, secret_keys);

    uint64_t ts = 0xDEADBEEF;
    pk = new Polynomial(agg->public_key(ts));

    std::cout << "Generating test data..." << std::endl;
    generate_all_test_data();
    std::cout << "Benchmark initialized." << std::endl;
}

TERSEBenchmark::~TERSEBenchmark() {
    delete agg_key;
    delete pk;
    delete agg;
}

void TERSEBenchmark::run_round(uint64_t round_idx,
                               double& sum_time_us,
                               double& lookup_time_us) {
    std::vector<Polynomial> ctexts;
    ctexts.reserve(config.num_devices);

    auto sum_start = get_time_us();

    // Encryption
    for(size_t i = 0; i < config.num_devices; i++) {
        Polynomial input(plain_parms);
        input.zero();
        input.at(0, 0) = all_user_inputs[round_idx][i] % plain_parms->moduli(0);

        double noise_time, enc_time;
        Polynomial ct = agg->enc(input, secret_keys[i], *pk,
                                true, noise_time, enc_time);
        ctexts.push_back(ct);
    }

    // Explicit aggregation (to match hybrid)
    Polynomial aggregated(ctext_parms);
    aggregated.zero();
    for(const auto& ct : ctexts) {
        aggregated += ct;
    }

    auto sum_end = get_time_us();
    sum_time_us = sum_end - sum_start;

    auto lookup_start = get_time_us();

    // Decryption only (single aggregated ciphertext)
    std::vector<Polynomial> agg_vec = {aggregated};
    double dec_time;
    uint64_t ts = 0xDEADBEEF;
    Polynomial result = agg->dec(*agg_key, agg_vec, ts, dec_time, 1);

    auto lookup_end = get_time_us();
    lookup_time_us = lookup_end - lookup_start;
}

void TERSEBenchmark::run_hybrid_native(uint64_t round_idx,
                                      HybridState& state,
                                      double& sum_time_us) {
    std::vector<Polynomial> ctexts;
    ctexts.reserve(config.num_devices);

    auto sum_start = get_time_us();

    // Encrypt all user inputs
    for(size_t i = 0; i < config.num_devices; i++) {
        Polynomial input(plain_parms);
        input.zero();
        input.at(0, 0) = all_user_inputs[round_idx][i] % plain_parms->moduli(0);

        double noise_time, enc_time;
        Polynomial ct = agg->enc(input, secret_keys[i], *pk,
                                true, noise_time, enc_time);
        ctexts.push_back(ct);
    }

    // Aggregate ciphertexts (outside SGX)
    Polynomial aggregated(ctext_parms);
    aggregated.zero();
    for(const auto& ct : ctexts) {
        aggregated += ct;
    }

    auto sum_end = get_time_us();
    sum_time_us = sum_end - sum_start;

    // Serialize aggregated ciphertext for SGX
    state.data_size = aggregated.size_in_bytes();
    state.serialized_agg_ctext = (uint64_t*)malloc(state.data_size);
    if (!state.serialized_agg_ctext) {
        throw std::runtime_error("Failed to allocate memory for serialization");
    }
    memcpy(state.serialized_agg_ctext, aggregated.buffer(), state.data_size);
}

void TERSEBenchmark::run_hybrid_sgx(const HybridState& state,
                                   double& lookup_time_us) {
    if (!state.serialized_agg_ctext) {
        throw std::runtime_error("Invalid hybrid state");
    }

    // Deserialize aggregated ciphertext
    Polynomial aggregated(ctext_parms, state.serialized_agg_ctext);

    auto lookup_start = get_time_us();

    // Key lookup and decryption inside SGX
    std::vector<Polynomial> ctexts = {aggregated};
    double dec_time;
    uint64_t ts = 0xDEADBEEF;
    Polynomial result = agg->dec(*agg_key, ctexts, ts, dec_time, 1);

    auto lookup_end = get_time_us();
    lookup_time_us = lookup_end - lookup_start;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <platform> <num_devices> <records_per_device>"
                  << std::endl;
        return 1;
    }

    std::string platform = argv[1];
    uint64_t num_devices = std::stoull(argv[2]);
    uint64_t records_per_device = std::stoull(argv[3]);

    BenchmarkConfig config;
    config.num_devices = num_devices;
    config.records_per_device = records_per_device;
    config.plain_bits = 32;
    config.scale = 0.5f;
    config.scheme = NS;
    config.platform = platform;

    TERSEBenchmark benchmark(config);

    std::string results_file = "terse_results.csv";
    bool file_exists = std::ifstream(results_file).good();

    std::ofstream outfile(results_file, std::ios::app);
    if (!file_exists) {
        outfile << "platform,num_devices,records_per_device,"
                << "avg_sum_per_round_us,avg_lookup_per_round_us,"
                << "avg_total_per_round_us" << std::endl;
    }

    if (platform == "native" || platform == "sgx") {
        double total_sum = 0.0;
        double total_lookup = 0.0;

        for(uint64_t round = 0; round < records_per_device; round++) {
            double sum_time, lookup_time;
            benchmark.run_round(round, sum_time, lookup_time);
            total_sum += sum_time;
            total_lookup += lookup_time;
        }

        double avg_sum = total_sum / records_per_device;
        double avg_lookup = total_lookup / records_per_device;
        double avg_total = avg_sum + avg_lookup;

        outfile << platform << "," << num_devices << ","
                << records_per_device << "," << avg_sum << ","
                << avg_lookup << "," << avg_total << std::endl;

        std::cout << "\n=== Results ===" << std::endl;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Avg sum: " << avg_sum << " us" << std::endl;
        std::cout << "Avg lookup: " << avg_lookup << " us" << std::endl;
        std::cout << "Avg total: " << avg_total << " us" << std::endl;

    } else if (platform == "hybrid_native") {
        double total_sum = 0.0;

        for(uint64_t round = 0; round < records_per_device; round++) {
            HybridState state;
            double sum_time;
            benchmark.run_hybrid_native(round, state, sum_time);
            total_sum += sum_time;
        }

        double avg_sum = total_sum / records_per_device;

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Sum: " << avg_sum << std::endl;

    } else if (platform == "hybrid_sgx") {
        double total_lookup = 0.0;

        for(uint64_t round = 0; round < records_per_device; round++) {
            HybridState state;
            Polynomial dummy_agg(benchmark.get_ctext_parms());
            dummy_agg.zero();

            state.data_size = dummy_agg.size_in_bytes();
            state.serialized_agg_ctext = (uint64_t*)malloc(state.data_size);
            memcpy(state.serialized_agg_ctext, dummy_agg.buffer(), state.data_size);

            double lookup_time;
            benchmark.run_hybrid_sgx(state, lookup_time);
            total_lookup += lookup_time;
        }

        double avg_lookup = total_lookup / records_per_device;

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Lookup: " << avg_lookup << std::endl;

    } else if (platform == "hybrid") {
        double total_sum = 0.0;
        double total_lookup = 0.0;

        for(uint64_t round = 0; round < records_per_device; round++) {
            HybridState state;
            double sum_time, lookup_time;

            benchmark.run_hybrid_native(round, state, sum_time);
            benchmark.run_hybrid_sgx(state, lookup_time);

            total_sum += sum_time;
            total_lookup += lookup_time;
        }

        double avg_sum = total_sum / records_per_device;
        double avg_lookup = total_lookup / records_per_device;
        double avg_total = avg_sum + avg_lookup;

        outfile << platform << "," << num_devices << ","
                << records_per_device << "," << avg_sum << ","
                << avg_lookup << "," << avg_total << std::endl;

        std::cout << "\n=== Hybrid Results ===" << std::endl;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Avg sum (native): " << avg_sum << " us" << std::endl;
        std::cout << "Avg lookup (SGX): " << avg_lookup << " us" << std::endl;
        std::cout << "Avg total: " << avg_total << " us" << std::endl;

    } else {
        std::cerr << "Unknown platform: " << platform << std::endl;
        return 1;
    }

    return 0;
}
