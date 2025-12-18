#!/bin/bash
set -e

DEVICE_COUNTS=(10000 100000 1000000)
RECORD_COUNTS=(15)

BENCHMARKS=(
    "aes:aes:native,sgx"
    "terse:terse:native,sgx,hybrid"
    "rsa:rsa:native,sgx"
    "ecc:ecc:native,sgx"
)

usage() {
    echo "Usage: $0 <command> [args...]"
    echo
    echo "Commands:"
    echo "  build <benchmark> <platform>           Build a specific benchmark"
    echo "  run <benchmark> <platform> <dev> <rec> Run a single benchmark"
    echo "  loop <benchmark> <platform>            Loop through all configs"
    echo "  loop-all <benchmark>                   Loop all platforms for benchmark"
    echo "  loop-everything                        Run all benchmarks, all platforms"
    echo "  aggregate                              Aggregate all results into one CSV"
    echo "  clean <benchmark> [all]                Clean artifacts"
    echo "  clean-all [all]                        Clean all benchmarks"
    echo
    echo "Benchmarks: aes, terse, rsa, ecc"
    echo "Platforms: native, sgx, sgx_only (terse), hybrid (terse)"
    exit 1
}

get_benchmark_info() {
    local bench=$1
    for entry in "${BENCHMARKS[@]}"; do
        IFS=':' read -r dir prefix platforms <<< "$entry"
        if [ "$dir" == "$bench" ]; then
            echo "$dir:$prefix:$platforms"
            return 0
        fi
    done
    echo "Error: Unknown benchmark '$bench'" >&2
    return 1
}

build_benchmark() {
    local bench=$1
    local platform=$2

    local info=$(get_benchmark_info "$bench")
    IFS=':' read -r dir prefix platforms <<< "$info"

    cd "$dir"

    local program="${prefix}_benchmark"
    local manifest_template="${prefix}_benchmark.manifest.template"
    local manifest="${prefix}_benchmark.manifest"
    local manifest_sgx="${prefix}_benchmark.manifest.sgx"

    echo "Building $bench for $platform..."

case $platform in
    native)
        if [ "$prefix" == "terse" ]; then
            g++ -std=c++17 -O3 -DNDEBUG -march=native -o "$program" \
                "${prefix}_benchmark.cpp" \
                -lntl -lgmp -lgmpxx -pthread
        elif [ "$prefix" == "aes" ] || [ "$prefix" == "rsa" ] || [ "$prefix" == "ecc" ]; then
            g++ -std=c++17 -O3 -DNDEBUG -march=native -o "$program" \
                "${prefix}_benchmark.cpp" \
                -lssl -lcrypto
        else
            g++ -std=c++17 -O3 -DNDEBUG -march=native -o "$program" \
                "${prefix}_benchmark.cpp"
        fi
        ;;
    sgx|sgx_only|hybrid)
        if [ "$prefix" == "terse" ]; then
            g++ -std=c++17 -O3 -DNDEBUG -march=native -o "$program" \
                "${prefix}_benchmark.cpp" \
                -lntl -lgmp -lgmpxx -pthread
        elif [ "$prefix" == "aes" ] || [ "$prefix" == "rsa" ] || [ "$prefix" == "ecc" ]; then
            g++ -std=c++17 -O3 -DNDEBUG -march=native -o "$program" \
                "${prefix}_benchmark.cpp" \
                -lssl -lcrypto
        else
            g++ -std=c++17 -O3 -DNDEBUG -march=native -o "$program" \
                "${prefix}_benchmark.cpp"
        fi
        gramine-manifest "$manifest_template" "$manifest"
        gramine-sgx-sign --manifest "$manifest" --output "$manifest_sgx"
        ;;
esac
    cd ..
    echo "Build complete: $bench ($platform)"
}

run_single() {
    local bench=$1
    local platform=$2
    local devices=$3
    local records=$4

    local info=$(get_benchmark_info "$bench")
    IFS=':' read -r dir prefix platforms <<< "$info"

    cd "$dir"

    local program="${prefix}_benchmark"
    local manifest_sgx="${prefix}_benchmark.manifest.sgx"

    if [ ! -f "$program" ]; then
        echo "Error: Executable not found. Build first." >&2
        cd ..
        exit 1
    fi

    if [[ "$platform" == "sgx" || "$platform" == "sgx_only" || "$platform" == "hybrid" ]] && [ ! -f "$manifest_sgx" ]; then
        echo "Error: SGX manifest not found. Build first." >&2
        cd ..
        exit 1
    fi

    echo "Running $bench on $platform: $devices devices, $records records..."

    case $platform in
        native)
            ./"$program" native "$devices" "$records"
            ;;
        sgx)
            gramine-sgx "$program" sgx "$devices" "$records"
            ;;
        sgx_only)
            gramine-sgx "$program" sgx_only "$devices" "$records"
            ;;
        hybrid)
            if [ "$prefix" != "terse" ]; then
                echo "Error: Hybrid mode only available for terse benchmarks" >&2
                cd ..
                exit 1
            fi
            run_hybrid "$program" "$devices" "$records"
            ;;
        *)
            echo "Error: Invalid platform '$platform'" >&2
            cd ..
            exit 1
            ;;
    esac

    cd ..
}

run_hybrid() {
    local program=$1
    local devices=$2
    local records=$3

    echo "--- Running HYBRID benchmark ---"

    echo "Step 1: Running native summation..."
    native_output=$(./"$program" hybrid_native "$devices" "$records" 2>&1)
    echo "$native_output"

    native_sum_time=$(echo "$native_output" | grep "Sum:" | awk '{print $2}')

    echo "Step 2: Running SGX lookup..."
    sgx_output=$(gramine-sgx "$program" hybrid_sgx "$devices" "$records" 2>&1)
    echo "$sgx_output"

    sgx_lookup_time=$(echo "$sgx_output" | grep "Lookup:" | awk '{print $2}')

    total_sum=$native_sum_time
    total_lookup=$sgx_lookup_time
    total_time=$(echo "scale=4; $total_sum + $total_lookup" | bc)

    echo ""
    echo "--- Hybrid Summary (avg per round) ---"
    echo "Sum:    $total_sum us"
    echo "Lookup: $total_lookup us"
    echo "Total:  $total_time us"

    local results_file="terse_results.csv"

    if [ ! -f "$results_file" ]; then
        echo "platform,num_devices,records_per_device,avg_sum_per_round_us,avg_lookup_per_round_us,avg_total_per_round_us" > "$results_file"
    fi
    echo "hybrid,$devices,$records,$total_sum,$total_lookup,$total_time" >> "$results_file"
}

loop_benchmark() {
    local bench=$1
    local platform=$2

    local info=$(get_benchmark_info "$bench")
    IFS=':' read -r dir prefix platforms <<< "$info"

    echo "====================================================="
    echo "Starting benchmark loop: $bench on $platform"
    echo "====================================================="

    build_benchmark "$bench" "$platform"

    for devices in "${DEVICE_COUNTS[@]}"; do
        for records in "${RECORD_COUNTS[@]}"; do
            echo
            echo "--- Config: $devices devices, $records records/device ---"
            run_single "$bench" "$platform" "$devices" "$records"
            echo "-------------------------------------------------------------"
        done
    done

    echo
    echo "====================================================="
    echo "Benchmark loop complete: $bench on $platform"
    echo "====================================================="
}

loop_all_platforms() {
    local bench=$1

    local info=$(get_benchmark_info "$bench")
    IFS=':' read -r dir prefix platforms <<< "$info"

    echo "=========================================="
    echo "Running all platforms for: $bench"
    echo "=========================================="

    IFS=',' read -ra platform_array <<< "$platforms"
    for platform in "${platform_array[@]}"; do
        echo
        echo "########## PLATFORM: $platform ##########"
        loop_benchmark "$bench" "$platform"
    done

    echo
    echo "=========================================="
    echo "All platforms complete for: $bench"
    echo "Results in $dir/"
    echo "=========================================="
}

loop_everything() {
    echo "############################################################"
    echo "# RUNNING ALL BENCHMARKS ON ALL PLATFORMS"
    echo "############################################################"
    echo

    for entry in "${BENCHMARKS[@]}"; do
        IFS=':' read -r dir prefix platforms <<< "$entry"
        echo
        echo "============================================================"
        echo "= BENCHMARK: $dir"
        echo "============================================================"
        loop_all_platforms "$dir"
    done

    echo
    echo "############################################################"
    echo "# ALL BENCHMARKS COMPLETE!"
    echo "############################################################"
    echo
    echo "Results locations:"
    for entry in "${BENCHMARKS[@]}"; do
        IFS=':' read -r dir prefix platforms <<< "$entry"
        echo "  - $dir/${prefix}_results.csv"
    done
}

aggregate_results() {
    local output_file="aggregated_results.csv"

    echo "Aggregating results into $output_file..."

    # Create header
    echo "benchmark,platform,num_devices,records_per_device,setup_us,avg_decrypt_sum_per_round_us,avg_final_per_round_us,avg_total_per_round_us" > "$output_file"

    # Process each benchmark
    for entry in "${BENCHMARKS[@]}"; do
        IFS=':' read -r dir prefix platforms <<< "$entry"
        local results_file="$dir/${prefix}_results.csv"

        if [ -f "$results_file" ]; then
            if [ "$prefix" == "terse" ]; then
                # Terse has different format: platform,devices,records,sum,lookup,total
                # Map to: benchmark,platform,devices,records,0,sum,lookup,total
                tail -n +2 "$results_file" | while IFS=, read -r platform devices records sum lookup total; do
                    echo "$prefix,$platform,$devices,$records,0,$sum,$lookup,$total" >> "$output_file"
                done
            else
                # Crypto benchmarks: platform,devices,records,setup,decrypt,final,total
                tail -n +2 "$results_file" | while IFS=, read -r platform devices records setup decrypt final total; do
                    echo "$prefix,$platform,$devices,$records,$setup,$decrypt,$final,$total" >> "$output_file"
                done
            fi
        else
            echo "Warning: $results_file not found, skipping..."
        fi
    done

    echo "Aggregation complete: $output_file"
    echo "Total rows: $(tail -n +2 "$output_file" | wc -l)"
}


clean_benchmark() {
    local bench=$1
    local clean_all=$2

    local info=$(get_benchmark_info "$bench")
    IFS=':' read -r dir prefix platforms <<< "$info"

    cd "$dir"

    local program="${prefix}_benchmark"
    local manifest="${prefix}_benchmark.manifest"
    local manifest_sgx="${prefix}_benchmark.manifest.sgx"

    if [ "$clean_all" == "all" ]; then
        echo "Cleaning all artifacts and results for $bench..."
        rm -f "$program" "$manifest" "$manifest_sgx"
        rm -f "${prefix}"*.csv
    else
        echo "Cleaning build artifacts for $bench..."
        rm -f "$program" "$manifest" "$manifest_sgx"
    fi

    cd ..
}

clean_all_benchmarks() {
    local clean_all=$1

    echo "Cleaning all benchmarks..."
    for entry in "${BENCHMARKS[@]}"; do
        IFS=':' read -r dir prefix platforms <<< "$entry"
        clean_benchmark "$dir" "$clean_all"
    done
    echo "Clean complete."
}

COMMAND=${1:-}
if [ -z "$COMMAND" ]; then usage; fi

case $COMMAND in
    build)
        if [ "$#" -ne 3 ]; then
            echo "Error: 'build' requires <benchmark> <platform>" >&2
            usage
        fi
        build_benchmark "$2" "$3"
        ;;
    run)
        if [ "$#" -ne 5 ]; then
            echo "Error: 'run' requires <benchmark> <platform> <devices> <records>" >&2
            usage
        fi
        run_single "$2" "$3" "$4" "$5"
        ;;
    loop)
        if [ "$#" -ne 3 ]; then
            echo "Error: 'loop' requires <benchmark> <platform>" >&2
            usage
        fi
        loop_benchmark "$2" "$3"
        ;;
    loop-all)
        if [ "$#" -ne 2 ]; then
            echo "Error: 'loop-all' requires <benchmark>" >&2
            usage
        fi
        loop_all_platforms "$2"
        ;;
    loop-everything)
        loop_everything
        ;;
    aggregate)
        aggregate_results
        ;;
    clean)
        if [ "$#" -lt 2 ]; then
            echo "Error: 'clean' requires <benchmark>" >&2
            usage
        fi
        clean_benchmark "$2" "${3:-}"
        ;;
    clean-all)
        clean_all_benchmarks "${2:-}"
        ;;
    *)
        echo "Error: Unknown command '$COMMAND'" >&2
        usage
        ;;
esac

echo "Done."
