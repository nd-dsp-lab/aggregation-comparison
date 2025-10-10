#!/bin/bash

set -e

# --- Configuration ---
PROGRAM_NAME="benchmark"
MANIFEST_TEMPLATE="benchmark.manifest.template"
MANIFEST="benchmark.manifest"
MANIFEST_SGX="benchmark.manifest.sgx"
DATA_FILE="encrypted_data.bin"
KEY_FILE="encrypted_keys.bin"
TERSE_DATA_FILE="terse_data.bin"
RESULTS_FILE="results.csv"

# --- Benchmark Loop Configuration ---
# Define the combinations of devices and records you want to test.
# Modify these arrays to control the benchmarks run by the 'loop-benchmark' command.
DEVICE_COUNTS=(1000 10000 100000)
RECORD_COUNTS=(1000)


# --- Helper Functions ---
usage() {
    echo "Usage: $0 <command> [args...]"
    echo
    echo "Commands:"
    echo "  all <platform> <devices> <records>  Build, generate, and run a single benchmark."
    echo "  build <platform>                    Build the benchmark executable."
    echo "  generate <platform> <devices> <records> Generate test data."
    echo "  run <platform> <devices> <records>  Run a benchmark."
    echo "  loop-benchmark <platform>           Loop through benchmarks defined in the script."
    echo "  clean [all]                         Remove build artifacts. 'all' also removes results."
    echo
    echo "Platforms:"
    echo "  native, sgx, terse-native, terse-sgx, terse-hybrid"
    exit 1
}

build_native() {
    echo "Building native binary..."
    g++ -std=c++17 -O3 -DNDEBUG -march=native -o "$PROGRAM_NAME" main.cpp -lssl -lcrypto
}

build_sgx() {
    echo "Building for SGX..."
    build_native
    gramine-manifest "$MANIFEST_TEMPLATE" "$MANIFEST"
    gramine-sgx-sign --manifest "$MANIFEST" --output "$MANIFEST_SGX"
}

# --- Action Functions ---
do_build() {
    local platform=$1
    case $platform in
        native|terse-native) build_native ;;
        sgx|terse-sgx|terse-hybrid) build_sgx ;;
        *) echo "Error: Invalid platform '$platform'." >&2; usage ;;
    esac
}

do_generate() {
    local platform=$1
    local devices=$2
    local records=$3
    if [ ! -f "$PROGRAM_NAME" ]; then echo "Error: Executable not found. Build first." >&2; exit 1; fi
    echo "--> Generating test data for $devices devices, $records records each..."
    case $platform in
        native|sgx) ./"$PROGRAM_NAME" generate "$devices" "$records" ;;
        terse-native|terse-sgx|terse-hybrid) ./"$PROGRAM_NAME" generate-terse "$devices" "$records" ;;
        *) echo "Error: Invalid platform '$platform' for generate." >&2; usage ;;
    esac
}

do_run() {
    local platform=$1
    local devices=$2
    local records=$3
    if [ ! -f "$PROGRAM_NAME" ]; then echo "Error: Executable not found. Build first." >&2; exit 1; fi
    if [[ "$platform" == "sgx" || "$platform" == "terse-sgx" || "$platform" == "terse-hybrid" ]] && [ ! -f "$MANIFEST_SGX" ]; then
        echo "Error: SGX manifest not found. Build first." >&2; exit 1;
    fi

    echo "--> Running benchmark for $devices devices, $records records each on platform: $platform"
    case $platform in
        native) ./"$PROGRAM_NAME" native "$devices" "$records" ;;
        sgx) gramine-sgx "$PROGRAM_NAME" sgx "$devices" "$records" ;;
        terse-native) ./"$PROGRAM_NAME" terse-native "$devices" "$records" ;;
        terse-sgx) gramine-sgx "$PROGRAM_NAME" terse-sgx "$devices" "$records" ;;
        terse-hybrid)
            echo "--- Measuring Part 1: Native Summation (Looped) ---"
            NATIVE_OUTPUT=$(./"$PROGRAM_NAME" terse-hybrid-native-sum-looped "$devices" "$records")
            NATIVE_TIME_US=$(echo "$NATIVE_OUTPUT" | grep 'duration_us:' | cut -d':' -f2)
            if [ -z "$NATIVE_TIME_US" ]; then echo "Error: Failed to get native duration."; exit 1; fi
            echo "Native part took: $NATIVE_TIME_US microseconds"

            echo "--- Measuring Part 2: SGX Lookup (Looped) ---"
            SGX_OUTPUT=$(gramine-sgx "$PROGRAM_NAME" terse-hybrid-sgx-lookup-looped "$devices" "$records")
            SGX_TIME_US=$(echo "$SGX_OUTPUT" | grep 'duration_us:' | cut -d':' -f2)
            if [ -z "$SGX_TIME_US" ]; then echo "Error: Failed to get SGX duration."; exit 1; fi
            echo "SGX part took: $SGX_TIME_US microseconds"

            TOTAL_TIME_US=$((NATIVE_TIME_US + SGX_TIME_US))
            TOTAL_RECORDS=$((devices * records))

            # Use bc for floating point calculations
            OPS_PER_SEC=$(echo "scale=2; if ($TOTAL_TIME_US > 0) $TOTAL_RECORDS * 1000000 / $TOTAL_TIME_US else 0" | bc)
            AVG_ROUND_TIME_US=$(echo "scale=2; if ($records > 0) $TOTAL_TIME_US / $records else 0" | bc)

            echo "----------------------------------------"
            echo "Total Time: $TOTAL_TIME_US microseconds"
            echo "Average round time: $AVG_ROUND_TIME_US microseconds"
            echo "Rate: $OPS_PER_SEC ops/sec"

            if [ ! -f "$RESULTS_FILE" ]; then
                echo "platform,num_devices,records_per_device,duration_us,ops_per_sec,avg_round_time_us" > "$RESULTS_FILE"
            fi
            echo "$platform,$devices,$records,$TOTAL_TIME_US,$OPS_PER_SEC,$AVG_ROUND_TIME_US" >> "$RESULTS_FILE"
            echo "Results appended to $RESULTS_FILE"
            ;;
        *) echo "Error: Invalid platform '$platform'." >&2; usage ;;
    esac
}

# --- Main Logic ---
COMMAND=${1:-}
PLATFORM=${2:-}
if [ -z "$COMMAND" ]; then usage; fi

case $COMMAND in
    all)
        if [ "$#" -ne 4 ]; then echo "Error: 'all' command requires <platform> <devices> <records>." >&2; usage; fi
        do_build "$2"
        do_generate "$2" "$3" "$4"
        do_run "$2" "$3" "$4"
        ;;
    build)
        if [ -z "$PLATFORM" ]; then echo "Error: 'build' command requires a platform." >&2; usage; fi
        do_build "$PLATFORM"
        ;;
    generate)
        if [ "$#" -ne 4 ]; then echo "Error: 'generate' command requires <platform> <devices> <records>." >&2; usage; fi
        do_generate "$2" "$3" "$4"
        ;;
    run)
        if [ "$#" -ne 4 ]; then echo "Error: 'run' command requires <platform> <devices> <records>." >&2; usage; fi
        do_run "$2" "$3" "$4"
        ;;
    clean)
        if [ "$2" == "all" ]; then
            echo "Cleaning all artifacts, data, and results..."
            rm -f "$PROGRAM_NAME" "$MANIFEST" "$MANIFEST_SGX" "$DATA_FILE" "$KEY_FILE" "$TERSE_DATA_FILE" "$RESULTS_FILE"
        else
            echo "Cleaning artifacts and data..."
            rm -f "$PROGRAM_NAME" "$MANIFEST" "$MANIFEST_SGX" "$DATA_FILE" "$KEY_FILE" "$TERSE_DATA_FILE"
        fi
        ;;
    loop-benchmark)
        if [ -z "$PLATFORM" ]; then
            echo "Error: 'loop-benchmark' command requires a platform." >&2
            usage
        fi
        echo "Starting full benchmark loop for platform: $PLATFORM"
        echo "====================================================="
        # Build once at the start
        do_build "$PLATFORM"

        # Loop through all combinations
        for devices in "${DEVICE_COUNTS[@]}"; do
            for records in "${RECORD_COUNTS[@]}"; do
                echo
                echo "--- BENCHMARKING: $devices devices, $records records/device ---"
                do_generate "$PLATFORM" "$devices" "$records"
                do_run "$PLATFORM" "$devices" "$records"
                echo "-------------------------------------------------------------"
            done
        done
        echo
        echo "====================================================="
        echo "Full benchmark loop finished. Results are in $RESULTS_FILE"
        ;;
    *)
        echo "Error: Unknown command '$COMMAND'." >&2; usage
        ;;
esac

echo "Done."
