# Benchmark Results

## Overview

Performance comparison of AES encryption and basic aggregation workloads across different execution environments and I/O methods. All measurements are in microseconds (μs) averaged per round over 1000 rounds.

**Platforms:**
- `native`: Standard Linux execution
- `sgx`: Intel SGX enclave via Gramine
- `hybrid`: Split execution (native aggregation + SGX lookup)

**I/O Methods:**
- `file`: Reading data from disk files
- `shm`: Reading from shared memory (untrusted memory region for SGX)

**Phases:**
- `avg_read_us`: One-time setup cost (loading keys/data)
- `avg_aggregation_us` / `avg_sum_us`: Decryption and summation per round
- `avg_write_us` / `avg_lookup_us`: Final aggregation or lookup table operations
- `avg_total_us`: Total time per round

## AES Benchmark

Encrypted data processing with per-device AES-128-CBC decryption and aggregation.

| benchmark | io_method | platform | num_devices | records_per_device | avg_read_us | avg_aggregation_us | avg_write_us | avg_total_us |
|-----------|-----------|----------|-------------|-------------------|-------------|-------------------|--------------|--------------|
| aes | file | native | 1000 | 1000 | 2004 | 603.457 | 0.001 | 603.458 |
| aes | file | native | 10000 | 1000 | 1847 | 6495.9 | 0.001 | 6495.9 |
| aes | file | native | 100000 | 1000 | 9323 | 72514.1 | 0.001 | 72514.1 |
| aes | file | sgx | 1000 | 1000 | 3133 | 812.91 | 0.002 | 812.912 |
| aes | file | sgx | 10000 | 1000 | 3169 | 6960.91 | 0.004 | 6960.92 |
| aes | file | sgx | 100000 | 1000 | 10746 | 81087.1 | 0.003 | 81087.1 |
| aes | shm | native | 1000 | 1000 | 1208 | 527.203 | 0.001 | 527.204 |
| aes | shm | native | 10000 | 1000 | 3108 | 6322.06 | 0.001 | 6322.06 |
| aes | shm | native | 100000 | 1000 | 9609 | 71435.6 | 0.001 | 71435.6 |
| aes | shm | sgx | 1000 | 1000 | 88 | 597.647 | 0.003 | 597.65 |
| aes | shm | sgx | 10000 | 1000 | 948 | 6959.18 | 0.002 | 6959.18 |
| aes | shm | sgx | 100000 | 1000 | 10916 | 69414.9 | 0.002 | 69414.9 |

## Terse Benchmark

Simple summation with lookup table operations (minimal encryption overhead).

| benchmark | io_method | platform | num_devices | records_per_device | avg_read_us | avg_sum_us | avg_lookup_us | avg_total_us |
|-----------|-----------|----------|-------------|-------------------|-------------|-----------|---------------|--------------|
| terse | file | native | 1000 | 1000 | 0 | 1.228 | 0.007 | 1.235 |
| terse | file | native | 10000 | 1000 | 0 | 75.39 | 0.011 | 75.401 |
| terse | file | native | 100000 | 1000 | 0 | 524.206 | 0.008 | 524.214 |
| terse | file | sgx | 1000 | 1000 | 0 | 0.659 | 0.005 | 0.664 |
| terse | file | sgx | 10000 | 1000 | 0 | 107.119 | 0.006 | 111.394 |
| terse | file | sgx | 100000 | 1000 | 0 | 1009.89 | 0.009 | 1015.16 |
| terse | file | hybrid | 1000 | 1000 | 0.1550 | 0.9970 | 0.0030 | 1.1550 |
| terse | file | hybrid | 10000 | 1000 | 0.1320 | 75.6000 | 0.0020 | 75.7340 |
| terse | file | hybrid | 100000 | 1000 | 0.1950 | 509.9500 | 0.0050 | 510.1500 |
| terse | shm | native | 1000 | 1000 | 0.001 | 1.01 | 0.005 | 1.016 |
| terse | shm | native | 10000 | 1000 | 0 | 75.603 | 0.012 | 75.616 |
| terse | shm | native | 100000 | 1000 | 0 | 512.523 | 0.008 | 512.532 |
| terse | shm | sgx | 1000 | 1000 | 0.001 | 11.457 | 0.008 | 11.467 |
| terse | shm | sgx | 10000 | 1000 | 0 | 137.811 | 0.008 | 138.736 |
| terse | shm | sgx | 100000 | 1000 | 0.001 | 1490.11 | 0.016 | 1502.08 |
| terse | shm | hybrid | 1000 | 1000 | 0.1160 | 1.1540 | 0.0130 | 1.2830 |
| terse | shm | hybrid | 10000 | 1000 | 0.1000 | 73.1460 | 0.0130 | 73.2590 |
| terse | shm | hybrid | 100000 | 1000 | 0.0960 | 513.7020 | 0.0130 | 513.8110 |