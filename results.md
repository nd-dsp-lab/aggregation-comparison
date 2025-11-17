# Benchmark Results

## Overview

Performance comparison of AES encryption and basic aggregation workloads using shared memory I/O across different execution environments. All measurements are in microseconds (μs) averaged over 50 rounds of aggregation.

**Platforms:**
- `native`: Standard Linux execution
- `sgx`: Intel SGX enclave via Gramine
- `hybrid`: Split execution (native aggregation + SGX lookup)

**Test Configuration:**
- I/O Method: Shared Memory (shm)
- Aggregation Rounds: 50 per test

## AES Benchmark (Shared Memory I/O)

Encrypted data processing with per-device AES-128-CBC decryption and aggregation.

| Platform | Num Devices | Setup (μs) | Avg Decrypt+Sum (μs) | Avg Final (μs) | Avg Total (μs) |
|----------|-------------|------------|---------------------|----------------|----------------|
| native | 10,000 | 2,240 | 5,146.26 | 0 | 5,146.26 |
| native | 1,000,000 | 83,909 | 894,071 | 0.02 | 894,071 |
| sgx | 10,000 | 780 | 6,048.52 | 0 | 6,048.52 |
| sgx | 1,000,000 | 83,366 | 687,335 | 0 | 687,335 |

## Terse Benchmark (Shared Memory I/O)

TERSE summation with lookup table operation.

| Platform | Num Devices | Avg Sum (μs) | Avg Lookup (μs) | Avg Total (μs) |
|----------|-------------|--------------|-----------------|----------------|
| native | 10,000 | 11.7 | 0 | 11.72 |
| native | 1,000,000 | 4,655.66 | 0 | 4,655.68 |
| native | 100,000,000 | 525,150 | 0.02 | 525,150 |
| sgx | 10,000 | 118.3 | 0.02 | 118.32 |
| sgx | 1,000,000 | 14,131.5 | 0.04 | 14,293.3 |
| sgx | 100,000,000 | 1,402,050 | 0.04 | 1,402,770 |
| hybrid | 10,000 | 11.84 | 0.22 | 12.06 |
| hybrid | 1,000,000 | 5,757.54 | 0.24 | 5,757.78 |
| hybrid | 100,000,000 | 539,268.6 | 0.18 | 539,268.78 |
