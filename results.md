# Benchmark Results

Measures average round aggregation time in μs over 1,000 rounds of aggregation.

## Platform Descriptions

- **Native:** AES decryption and summation of each cipher text run without SGX.

- **SGX:** AES decryption and summation of each cipher text run inside SGX.

- **TERSE-Native:** TERSE run entirely without SGX.

- **TERSE-SGX:** TERSE run entirely in SGX.

- **TERSE-Hybrid (ours):** Most summations run outside of SGX. Lookup and final addition run inside SGX.

## Results

| Number of Devices | Native    | SGX       | TERSE-Native | TERSE-SGX | TERSE-Hybrid |
| ----------------- | --------- | --------- | ------------ | --------- | ------------ |
| 1,000             | 795.18    | 607.61    | 1.44         | 26.68     | 0.94         |
| 10,000            | 5,002.92  | 5,675.67  | 20.11        | 233.49    | 76.09        |
| 100,000           | 49,172.41 | 57,591.84 | 224.41       | 2,103.93  | 525.86       |

## Notes

- Results at 1,000 devices may exhibit timing artifacts due to measurement precision at microsecond scales. Focus on 10,000+ device results for reliable comparisons.
