#!/bin/bash

set -e

# Configuration
DEVICE_COUNTS=(1000 10000 100000)
RECORD_COUNTS=(1000)

BENCHMARKS=(
    "aes:aes:native,sgx"
    "aes-shm:aes:native,sgx"
    "terse:terse:native,sgx,hybrid"
    "terse-shm:terse:native,sgx,hybrid"
)

usage() {
    echo "Usage: $0 <command> [args...]"
    echo
    echo "Commands:"
    echo "  build <benchmark> <platform>           Build a specific benchmark"
    echo "  generate <benchmark> <devices> <recs>  Generate test data"
    echo "  run <benchmark> <platform> <dev> <rec> Run a single benchmark"
    echo "  loop <benchmark> <platform>            Loop through all configs"
    echo "  loop-all <benchmark>                   Loop all platforms for benchmark"
    echo "  loop-everything                        Run all benchmarks, all platforms"
    echo "  clean <benchmark> [all]                Clean artifacts"
    echo "  clean-all [all]                        Clean all benchmarks"
    echo "  aggregate                              Combine all results into one CSV"
    echo
    echo "Benchmarks: aes, aes-shm, terse, terse-shm"
    echo "Platforms: native, sgx, hybrid (terse only)"
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
            if [ "$prefix" == "aes" ]; then
                g++ -std=c++17 -O3 -DNDEBUG -march=native -o "$program" "${prefix}_benchmark.cpp" -lssl -lcrypto
            else
                g++ -std=c++17 -O3 -DNDEBUG -march=native -o "$program" "${prefix}_benchmark.cpp"
            fi
            ;;
        sgx|hybrid)
            if [ "$prefix" == "aes" ]; then
                g++ -std=c++17 -O3 -DNDEBUG -march=native -o "$program" "${prefix}_benchmark.cpp" -lssl -lcrypto
            else
                g++ -std=c++17 -O3 -DNDEBUG -march=native -o "$program" "${prefix}_benchmark.cpp"
            fi
            gramine-manifest "$manifest_template" "$manifest"
            gramine-sgx-sign --manifest "$manifest" --output "$manifest_sgx"
            ;;
        *)
            echo "Error: Invalid platform '$platform'" >&2
            cd ..
            exit 1
            ;;
    esac

    cd ..
    echo "Build complete: $bench ($platform)"
}

generate_data() {
    local bench=$1
    local devices=$2
    local records=$3

    local info=$(get_benchmark_info "$bench")
    IFS=':' read -r dir prefix platforms <<< "$info"

    cd "$dir"

    local program="${prefix}_benchmark"

    if [ ! -f "$program" ]; then
        echo "Error: Executable not found. Build first." >&2
        cd ..
        exit 1
    fi

    echo "Generating data for $bench: $devices devices, $records records..."
    ./"$program" generate "$devices" "$records"

    cd ..
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

    if [[ "$platform" == "sgx" || "$platform" == "hybrid" ]] && [ ! -f "$manifest_sgx" ]; then
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

    native_read_time=$(echo "$native_output" | grep "Read:" | awk '{print $2}')
    native_sum_time=$(echo "$native_output" | grep "Sum:" | awk '{print $2}')

    echo "Step 2: Running SGX lookup..."
    sgx_output=$(gramine-sgx "$program" hybrid_sgx "$devices" "$records" 2>&1)
    echo "$sgx_output"

    sgx_read_time=$(echo "$sgx_output" | grep "Read:" | awk '{print $2}')
    sgx_lookup_time=$(echo "$sgx_output" | grep "Lookup:" | awk '{print $2}')

    total_read=$(echo "scale=4; $native_read_time + $sgx_read_time" | bc)
    total_sum=$native_sum_time
    total_lookup=$sgx_lookup_time
    total_time=$(echo "scale=4; $total_read + $total_sum + $total_lookup" | bc)

    echo ""
    echo "--- Hybrid Summary (avg per round) ---"
    echo "Read:   $total_read us"
    echo "Sum:    $total_sum us"
    echo "Lookup: $total_lookup us"
    echo "Total:  $total_time us"

    local results_file="terse_results.csv"
    local breakdown_file="terse_hybrid_breakdown.csv"

    if [ ! -f "$results_file" ]; then
        echo "platform,num_devices,records_per_device,avg_read_per_round_us,avg_sum_per_round_us,avg_lookup_per_round_us,avg_total_per_round_us" > "$results_file"
    fi
    echo "hybrid,$devices,$records,$total_read,$total_sum,$total_lookup,$total_time" >> "$results_file"

    if [ ! -f "$breakdown_file" ]; then
        echo "num_devices,records_per_device,native_read_us,native_sum_us,sgx_read_us,sgx_lookup_us,total_us" > "$breakdown_file"
    fi
    echo "$devices,$records,$native_read_time,$native_sum_time,$sgx_read_time,$sgx_lookup_time,$total_time" >> "$breakdown_file"
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
            generate_data "$bench" "$devices" "$records"
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
        if [[ "$platforms" == *"hybrid"* ]]; then
            echo "  - $dir/${prefix}_hybrid_breakdown.csv"
        fi
    done
}

aggregate_results() {
    echo "Creating separate aggregated results..."

    # AES results
    local aes_file="all_aes_results.csv"
    echo "benchmark,io_method,platform,num_devices,records_per_device,avg_read_us,avg_aggregation_us,avg_write_us,avg_total_us" > "$aes_file"

    for dir in aes aes-shm; do
        # Results file is always named with base prefix (aes_results.csv in both dirs)
        results_file="$dir/aes_results.csv"

        if [ -f "$results_file" ]; then
            echo "  Adding $results_file to $aes_file..."
            io_method=$([ "$dir" == "aes-shm" ] && echo "shm" || echo "file")
            tail -n +2 "$results_file" | sed "s/^/aes,$io_method,/" >> "$aes_file"
        else
            echo "  Warning: $results_file not found"
        fi
    done

    # Terse results
    local terse_file="all_terse_results.csv"
    echo "benchmark,io_method,platform,num_devices,records_per_device,avg_read_us,avg_sum_us,avg_lookup_us,avg_total_us" > "$terse_file"

    for dir in terse terse-shm; do
        # Results file is always named with base prefix (terse_results.csv in both dirs)
        results_file="$dir/terse_results.csv"

        if [ -f "$results_file" ]; then
            echo "  Adding $results_file to $terse_file..."
            io_method=$([ "$dir" == "terse-shm" ] && echo "shm" || echo "file")
            tail -n +2 "$results_file" | sed "s/^/terse,$io_method,/" >> "$terse_file"
        else
            echo "  Warning: $results_file not found"
        fi
    done

    echo ""
    echo "Aggregation complete:"
    echo "  AES:   $aes_file ($(tail -n +2 "$aes_file" | wc -l) rows)"
    echo "  Terse: $terse_file ($(tail -n +2 "$terse_file" | wc -l) rows)"
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
        rm -f "${prefix}"*.bin "${prefix}"*.csv intermediate.bin
    else
        echo "Cleaning build artifacts for $bench..."
        rm -f "$program" "$manifest" "$manifest_sgx"
        rm -f "${prefix}"*.bin intermediate.bin
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

# Main command dispatcher
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
    generate)
        if [ "$#" -ne 4 ]; then
            echo "Error: 'generate' requires <benchmark> <devices> <records>" >&2
            usage
        fi
        generate_data "$2" "$3" "$4"
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
    aggregate)
        aggregate_results
        ;;
    *)
        echo "Error: Unknown command '$COMMAND'" >&2
        usage
        ;;
esac

echo "Done."
