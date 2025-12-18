#ifndef TERSE_BENCHMARK_H
#define TERSE_BENCHMARK_H

#include "Aggregator_RNS.h"
#include "DiscreteLaplacian.h"
#include "Polynomial.h"
#include "scheme.h"
#include <vector>
#include <string>
#include <chrono>

struct BenchmarkConfig {
    uint64_t num_devices;
    uint64_t records_per_device;
    unsigned int plain_bits;
    float scale;
    Scheme scheme;
    std::string platform;
};

struct HybridState {
    uint64_t* serialized_agg_ctext;
    size_t data_size;

    HybridState() : serialized_agg_ctext(nullptr), data_size(0) {}

    ~HybridState() {
        if (serialized_agg_ctext) {
            free(serialized_agg_ctext);
            serialized_agg_ctext = nullptr;
        }
    }
};

class TERSEBenchmark {
private:
    Aggregator_RNS* agg;
    Parameters* ctext_parms;
    Parameters* plain_parms;
    Polynomial* agg_key;
    std::vector<Polynomial> secret_keys;
    Polynomial* pk;
    BenchmarkConfig config;

    // Pre-generated test data for fair comparison
    std::vector<std::vector<uint64_t>> all_user_inputs;

    void generate_all_test_data();

public:
    TERSEBenchmark(const BenchmarkConfig& config);
    ~TERSEBenchmark();

    uint64_t get_time_us();

    Parameters* get_ctext_parms() const { return ctext_parms; }

    void run_round(uint64_t round_idx,
                   double& sum_time_us,
                   double& lookup_time_us);

    void run_hybrid_native(uint64_t round_idx,
                          HybridState& state,
                          double& sum_time_us);

    void run_hybrid_sgx(const HybridState& state,
                       double& lookup_time_us);
};

#endif
