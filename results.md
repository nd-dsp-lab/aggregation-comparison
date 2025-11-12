# Benchmark Results

## Overview

Performance comparison of AES encryption and basic aggregation workloads across different execution environments and I/O methods. All measurements are in microseconds (μs) averaged per round over 1000 rounds.

**Platforms:**
- `native`: Standard Linux execution
- `sgx`: Intel SGX enclave via Gramine
- `hybrid`: Split execution (native aggregation + SGX lookup)

**Phases:**
- `avg_read_us`: One-time setup cost (loading keys/data)
- `avg_aggregation_us` / `avg_sum_us`: Decryption and summation per round
- `avg_write_us` / `avg_lookup_us`: Final aggregation or lookup table operations
- `avg_total_us`: Total time per round

**Note:** All tests use 1000 records per device.

## AES Benchmark (File I/O)

Encrypted data processing with per-device AES-128-CBC decryption and aggregation.

| platform | num_devices | avg_read_us | avg_aggregation_us | avg_write_us | avg_total_us |
|----------|-------------|-------------|-------------------|--------------|--------------|
| native | 1000 | 2004 | 603.457 | 0.001 | 603.458 |
| native | 10000 | 1847 | 6495.9 | 0.001 | 6495.9 |
| native | 100000 | 9323 | 72514.1 | 0.001 | 72514.1 |
| sgx | 1000 | 3133 | 812.91 | 0.002 | 812.912 |
| sgx | 10000 | 3169 | 6960.91 | 0.004 | 6960.92 |
| sgx | 100000 | 10746 | 81087.1 | 0.003 | 81087.1 |

## Terse Benchmark (File I/O)

Simple summation with lookup table operations (minimal encryption overhead).

| platform | num_devices | avg_read_us | avg_sum_us | avg_lookup_us | avg_total_us |
|----------|-------------|-------------|-----------|---------------|--------------|
| native | 1000 | 0 | 1.228 | 0.007 | 1.235 |
| native | 10000 | 0 | 75.39 | 0.011 | 75.401 |
| native | 100000 | 0 | 524.206 | 0.008 | 524.214 |
| sgx | 1000 | 0 | 0.659 | 0.005 | 0.664 |
| sgx | 10000 | 0 | 107.119 | 0.006 | 111.394 |
| sgx | 100000 | 0 | 1009.89 | 0.009 | 1015.16 |
| hybrid | 1000 | 0.1550 | 0.9970 | 0.0030 | 1.1550 |
| hybrid | 10000 | 0.1320 | 75.6000 | 0.0020 | 75.7340 |
| hybrid | 100000 | 0.1950 | 509.9500 | 0.0050 | 510.1500 |

## AES Benchmark (Shared Memory I/O)

| platform | num_devices | avg_read_us | avg_aggregation_us | avg_write_us | avg_total_us |
|----------|-------------|-------------|-------------------|--------------|--------------|
| native | 1000 | 1208 | 527.203 | 0.001 | 527.204 |
| native | 10000 | 3108 | 6322.06 | 0.001 | 6322.06 |
| native | 100000 | 9609 | 71435.6 | 0.001 | 71435.6 |
| sgx | 1000 | 88 | 597.647 | 0.003 | 597.65 |
| sgx | 10000 | 948 | 6959.18 | 0.002 | 6959.18 |
| sgx | 100000 | 10916 | 69414.9 | 0.002 | 69414.9 |

## Terse Benchmark (Shared Memory I/O)

| platform | num_devices | avg_read_us | avg_sum_us | avg_lookup_us | avg_total_us |
|----------|-------------|-------------|-----------|---------------|--------------|
| native | 1000 | 0.001 | 1.01 | 0.005 | 1.016 |
| native | 10000 | 0 | 75.603 | 0.012 | 75.616 |
| native | 100000 | 0 | 512.523 | 0.008 | 512.532 |
| sgx | 1000 | 0.001 | 11.457 | 0.008 | 11.467 |
| sgx | 10000 | 0 | 137.811 | 0.008 | 138.736 |
| sgx | 100000 | 0.001 | 1490.11 | 0.016 | 1502.08 |
| hybrid | 1000 | 0.1160 | 1.1540 | 0.0130 | 1.2830 |
| hybrid | 10000 | 0.1000 | 73.1460 | 0.0130 | 73.2590 |
| hybrid | 100000 | 0.0960 | 513.7020 | 0.0130 | 513.8110 |
