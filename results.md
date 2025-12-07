# Benchmark Results

## Overview

Performance comparison of secure aggregation protocols in native and SGX. All measurements are in microseconds (μs) averaged per round.

**Platforms:**
- `native`: Standard Linux execution
- `sgx` and `sgx_only`: Intel SGX enclave via Gramine
- `hybrid`: Split execution - native aggregation + SGX lookup (TERSE)

**Test Configuration:**
- Records per Device: 15 (except RSA tests)
- Device Counts: 100,000 / 1,000,000 / 10,000,000

## AES Benchmark

AES-128-CBC encryption with per-device key decryption and aggregation.

| Platform | Num Devices | Setup (μs) | Avg Decrypt+Sum (μs) | Avg Final (μs) | Avg Total (μs) |
|----------|-------------|------------|---------------------|----------------|----------------|
| native | 100,000 | 7,094 | 63,187.8 | 0 | 63,187.8 |
| native | 1,000,000 | 74,542 | 848,538 | 0 | 848,538 |
| native | 10,000,000 | 795,882 | 10,354,700 | 0 | 10,354,700 |
| sgx | 100,000 | 8,435 | 68,355.9 | 0.067 | 68,356 |
| sgx | 1,000,000 | 100,354 | 926,692 | 0 | 926,692 |
| sgx | 10,000,000 | 818,230 | 12,314,200 | 0 | 12,314,200 |

## ECC Benchmark

Elliptic Curve Cryptography (ECIES) with P-256 curve encryption/decryption.

| Platform | Num Devices | Setup (μs) | Avg Decrypt+Sum (μs) | Avg Final (μs) | Avg Total (μs) |
|----------|-------------|------------|---------------------|----------------|----------------|
| native | 100,000 | 0 | 18,102,300 | 0 | 18,102,300 |
| native | 1,000,000 | 0 | 182,941,000 | 0.067 | 182,941,000 |
| native | 10,000,000 | 0 | 1,807,080,000 | 0 | 1,807,080,000 |
| sgx | 100,000 | 0 | 18,565,300 | 0.067 | 18,565,300 |
| sgx | 1,000,000 | 0 | 185,900,000 | 0.067 | 185,900,000 |
| sgx | 10,000,000 | 0 | 1,863,010,000 | 0 | 1,863,010,000 |


## RSA Benchmark

RSA-3072 encryption with OAEP padding. (Only 3 rounds)

| Platform | Num Devices | Setup (μs) | Avg Decrypt+Sum (μs) | Avg Final (μs) | Avg Total (μs) |
|----------|-------------|------------|---------------------|----------------|----------------|
| native | 100,000 | 0 | 68,852,000 | 0 | 68,852,000 |
| native | 1,000,000 | 0 | 688,378,000 | 0 | 688,378,000 |
| native | 10,000,000 | 1 | 6,880,870,000 | 0 | 6,880,870,000 |
| sgx | 100,000 | 0 | 70,405,300 | 0 | 70,405,300 |
| sgx | 1,000,000 | 0 | 703,962,000 | 0 | 703,962,000 |
| sgx | 10,000,000 | 0 | 7,038,560,000 | 0 | 7,038,560,000 |

## TERSE Benchmark

TERSE summation with lookup table.

| Platform | Num Devices | Avg Sum (μs) | Avg Lookup (μs) | Avg Total (μs) |
|----------|-------------|--------------|-----------------|----------------|
| native | 100,000 | 539.7 | 0.067 | 539.7 |
| native | 1,000,000 | 6,217.2 | 0 | 6,217.2 |
| native | 10,000,000 | 51,063.6 | 0.067 | 51,063.7 |
| native | 100,000,000 | 469,355 | 0.067 | 469,355 |
| sgx_only | 100,000 | 238.1 | 0 | 238.2 |
| sgx_only | 1,000,000 | 4,337.6 | 0.067 | 4,343.2 |
| sgx_only | 10,000,000 | 50,062.1 | 0.067 | 50,923.5 |
| sgx_only | 100,000,000 | 553,761 | 0.2 | 554,880 |
| hybrid | 100,000 | 397.3 | 0 | 397.3 |
| hybrid | 1,000,000 | 5,188.3 | 0 | 5,188.3 |
| hybrid | 10,000,000 | 46,926.0 | 0 | 46,926.0 |
| hybrid | 100,000,000 | 467,504.3 | 0.067 | 467,504.3 |