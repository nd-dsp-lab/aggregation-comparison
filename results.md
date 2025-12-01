# Benchmark Results

## Overview

Performance comparison of secure aggregation protocols in native and SGX. All measurements are in microseconds (μs) averaged per round.

**Platforms:**
- `native`: Standard Linux execution
- `sgx` and `sgx_only`: Intel SGX enclave via Gramine
- `hybrid`: Split execution - native aggregation + SGX lookup (TERSE)

**Test Configuration:**
- Records per Device: **15** (all tests)
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

RSA-2048 encryption with OAEP padding.

| Platform | Num Devices | Setup (μs) | Avg Decrypt+Sum (μs) | Avg Final (μs) | Avg Total (μs) |
|----------|-------------|------------|---------------------|----------------|----------------|
| native | 100,000 | 0 | 26,058,200 | 0 | 26,058,200 |
| native | 1,000,000 | 1 | 259,903,000 | 0 | 259,903,000 |
| native | 10,000,000 | 0 | 2,599,660,000 | 0 | 2,599,660,000 |
| sgx | 100,000 | 0 | 26,844,800 | 0.067 | 26,844,800 |
| sgx | 1,000,000 | 1 | 271,563,000 | 0 | 271,563,000 |
| sgx | 10,000,000 | 1 | 2,705,040,000 | 0 | 2,705,040,000 |

## TERSE Benchmark

TERSE summation with lookup table operation (non-cryptographic baseline).

| Platform | Num Devices | Avg Sum (μs) | Avg Lookup (μs) | Avg Total (μs) |
|----------|-------------|--------------|-----------------|----------------|
| native | 100,000 | 396.3 | 0.067 | 396.4 |
| native | 1,000,000 | 4,155.6 | 0 | 4,155.6 |
| native | 10,000,000 | 45,083.6 | 0.067 | 45,083.7 |
| sgx_only | 100,000 | 252.7 | 0.067 | 252.9 |
| sgx_only | 1,000,000 | 4,955.3 | 0 | 5,186.5 |
| sgx_only | 10,000,000 | 49,085.9 | 0.067 | 49,570.9 |
| hybrid | 100,000 | 227.0 | 0 | 227.0 |
| hybrid | 1,000,000 | 5,609.2 | 0 | 5,609.2 |
| hybrid | 10,000,000 | 45,335.1 | 0.067 | 45,335.2 |