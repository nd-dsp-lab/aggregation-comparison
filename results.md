# Benchmark Results

## Overview

Performance comparison of AES encryption and basic aggregation workloads across different execution environments using shared memory I/O. All measurements are in microseconds (μs) averaged per round over 1000 rounds.

**Platforms:**
- `native`: Standard Linux execution
- `sgx`: Intel SGX enclave via Gramine
- `hybrid`: Split execution (native aggregation + SGX lookup)

**Phases:**
- Setup: One-time initialization cost (loading keys/data)
- Per-round operations: Decryption, summation, and lookup operations
- Total: Complete per-round execution time

## AES Benchmark

Encrypted data processing with per-device AES-128-CBC decryption and aggregation using shared memory.

| Platform | Devices | Records/Device | Setup (μs) | Decrypt+Sum (μs) | Final (μs) | Total (μs) |
|----------|---------|----------------|------------|------------------|------------|------------|
| native | 1,000 | 1,000 | 1,208 | 527.203 | 0.001 | 527.204 |
| native | 10,000 | 1,000 | 3,108 | 6,322.06 | 0.001 | 6,322.06 |
| native | 10,000 | 50 | 3,160 | 7,844.78 | 0 | 7,844.78 |
| native | 100,000 | 1,000 | 9,609 | 71,435.6 | 0.001 | 71,435.6 |
| native | 1,000,000 | 50 | 83,304 | 892,824 | 0.02 | 892,824 |
| sgx | 1,000 | 1,000 | 88 | 597.647 | 0.003 | 597.65 |
| sgx | 10,000 | 1,000 | 948 | 6,959.18 | 0.002 | 6,959.18 |
| sgx | 10,000 | 50 | 819 | 5,863.86 | 0 | 5,863.86 |
| sgx | 100,000 | 1,000 | 10,916 | 69,414.9 | 0.002 | 69,414.9 |
| sgx | 1,000,000 | 50 | 99,074 | 709,714 | 0.02 | 709,714 |

## Terse Benchmark

Simple summation with lookup table operations using shared memory (minimal encryption overhead).

| Platform | Devices | Records/Device | Read (μs) | Sum (μs) | Lookup (μs) | Total (μs) |
|----------|---------|----------------|-----------|----------|-------------|------------|
| native | 1,000 | 1,000 | 0.001 | 1.01 | 0.005 | 1.016 |
| native | 10,000 | 1,000 | 0 | 75.603 | 0.012 | 75.616 |
| native | 10,000 | 50 | 0 | 12.48 | 0.02 | 12.5 |
| native | 100,000 | 1,000 | 0 | 512.523 | 0.008 | 512.532 |
| native | 1,000,000 | 50 | 0 | 4,655.02 | 0 | 4,655.04 |
| native | 100,000,000 | 50 | 0 | 524,403 | 0.02 | 524,403 |
| sgx | 1,000 | 1,000 | 0.001 | 11.457 | 0.008 | 11.467 |
| sgx | 10,000 | 1,000 | 0 | 137.811 | 0.008 | 138.736 |
| sgx | 10,000 | 50 | 0 | 120.18 | 0.02 | 120.2 |
| sgx | 100,000 | 1,000 | 0.001 | 1,490.11 | 0.016 | 1,502.08 |
| sgx | 1,000,000 | 50 | 0 | 13,694.1 | 0.04 | 13,709.4 |
| sgx | 100,000,000 | 50 | 0 | 1,417,580 | 0.02 | 1,417,590 |
| hybrid | 1,000 | 1,000 | 0.1160 | 1.1540 | 0.0130 | 1.2830 |
| hybrid | 10,000 | 1,000 | 0.1000 | 73.1460 | 0.0130 | 73.2590 |
| hybrid | 10,000 | 50 | 2.7600 | 11.6800 | 0.3000 | 14.7400 |
| hybrid | 100,000 | 1,000 | 0.0960 | 513.7020 | 0.0130 | 513.8110 |
| hybrid | 1,000,000 | 50 | 2.0600 | 4,618.6400 | 0.1800 | 4,620.8800 |
| hybrid | 100,000,000 | 50 | 2.2600 | 525,942.8600 | 0.1800 | 525,945.3000 |
