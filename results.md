Measures average round aggregation time in μs over 1,000 rounds of aggregation.

- **Native:** AES decryption and summation of each cipher text run without SGX.

- **SGX:** AES decryption and summation of each cipher text run inside SGX.

- **TERSE-Native:** TERSE run entirely without SGX.

- **TERSE-SGX:** TERSE run entirely in SGX.

- **TERSE-Hybrid (ours):** Most summations run outside of SGX. Lookup and final addition run inside SGX.

| Number of Devices | Native    | SGX       | Terse-Native | Terse-SGX | Terse-Hybrid |
| ----------------- | --------- | --------- | ------------ | --------- | ------------ |
| 1,000             | 630.07    | 578.57    | 1.04         | 25.27     | 1.90         |
| 10,000            | 4,792.59  | 5,773.88  | 13.70        | 253.41    | 15.46        |
| 100,000           | 47,380.04 | 57,990.95 | 194.53       | 1,985.00  | 144.76       |