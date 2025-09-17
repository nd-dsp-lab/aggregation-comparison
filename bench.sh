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
DEVICE_COUNTS=(10 100 1000 10000)
RECORD_COUNTS=(1 10 100 1000 10000)


# --- Helper Functions ---
usage() {
    echo "Usage: $0 <command> [args...]"
    echo
    echo "Commands:"
    echo "  all <platform> <devices> <records>  Build, generate, and run a single benchmark."
    echo "  build <platform>                    Build the benchmark executable."
    echo "  generate <platform> <devices> <records> Generate test data for a specific size."
    echo "  run <platform> <devices> <records>  Run a benchmark with pre-generated data."
    echo "  loop-benchmark <platform>           Build once, then loop through generating and running"
    echo "                                      benchmarks for all combinations in DEVICE_COUNTS"
    echo "                                      and RECORD_COUNTS defined in the script."
    echo "  clean                               Remove build artifacts and generated data."
    echo
    echo "Platforms:"
    echo "  native, sgx, terse-sgx"
    exit 1
}

build_native() {
    echo "Building native binary..."
    GRAMINE_API_FLAGS=$(pkg-config --cflags --libs gramine-api 2>/dev/null || echo "")
    g++ -std=c++17 -O3 -o "$PROGRAM_NAME" main.cpp -lssl -lcrypto $GRAMINE_API_FLAGS
}

build_sgx() {
    echo "Building for SGX..."
    build_native # SGX build requires the native binary first

    if [ ! -f "$MANIFEST_TEMPLATE" ]; then
        echo "Error: Manifest template '$MANIFEST_TEMPLATE' not found."
        exit 1
    fi

    echo "Generating manifest..."
    gramine-manifest "$MANIFEST_TEMPLATE" "$MANIFEST"

    echo "Signing enclave..."
    gramine-sgx-sign --manifest "$MANIFEST" --output "$MANIFEST_SGX"
}

# --- Action Functions ---
do_build() {
    local platform=$1
    case $platform in
        native) build_native ;;
        sgx|terse-sgx) build_sgx ;;
        *) echo "Error: Invalid platform '$platform'." >&2; usage ;;
    esac
}

do_generate() {
    local platform=$1
    local devices=$2
    local records=$3

    if [ ! -f "$PROGRAM_NAME" ]; then
        echo "Error: Executable not found. Please run 'build' first." >&2
        exit 1
    fi
    echo "--> Generating test data for $devices devices, $records records each..."
    # Generation can always run natively
    case $platform in
        native|sgx) ./"$PROGRAM_NAME" generate "$devices" "$records" ;;
        terse-sgx) ./"$PROGRAM_NAME" generate-terse "$devices" "$records" ;;
        *) echo "Error: Invalid platform '$platform' for generate." >&2; usage ;;
    esac
}

do_run() {
    local platform=$1
    local devices=$2
    local records=$3

    if [ ! -f "$PROGRAM_NAME" ]; then
        echo "Error: Executable not found. Please run 'build' first." >&2
        exit 1
    fi
    if [[ "$platform" == "sgx" || "$platform" == "terse-sgx" ]] && [ ! -f "$MANIFEST_SGX" ]; then
        echo "Error: Executable/manifest not found. Please run 'build' first." >&2
        exit 1
    fi

    # Check for correct data files
    if [[ "$platform" == "native" || "$platform" == "sgx" ]] && [[ ! -f "$DATA_FILE" || ! -f "$KEY_FILE" ]]; then
        echo "Error: Data files not found. Please run 'generate $platform <dev> <rec>' first." >&2
        exit 1
    elif [[ "$platform" == "terse-sgx" ]] && [[ ! -f "$TERSE_DATA_FILE" ]]; then
        echo "Error: Data file ($TERSE_DATA_FILE) not found. Please run 'generate $platform <dev> <rec>' first." >&2
        exit 1
    fi

    echo "--> Running benchmark for $devices devices, $records records each on platform: $platform"
    case $platform in
        native) ./"$PROGRAM_NAME" native "$devices" "$records" ;;
        sgx) gramine-sgx "$PROGRAM_NAME" sgx "$devices" "$records" ;;
        terse-sgx) gramine-sgx "$PROGRAM_NAME" terse-sgx "$devices" "$records" ;;
        *) echo "Error: Invalid platform '$platform'." >&2; usage ;;
    esac
}

# --- Main Logic ---
COMMAND=${1:-}
PLATFORM=${2:-}

if [ -z "$COMMAND" ]; then
    usage
fi

case $COMMAND in
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

    all)
        if [ "$#" -ne 4 ]; then
            echo "Error: 'all' command requires <platform> <devices> <records>." >&2
            usage
        fi
        do_build "$2"
        do_generate "$2" "$3" "$4"
        do_run "$2" "$3" "$4"
        ;;

    build)
        if [ -z "$PLATFORM" ]; then
            echo "Error: 'build' command requires a platform." >&2
            usage
        fi
        do_build "$PLATFORM"
        ;;

    generate)
        if [ "$#" -ne 4 ]; then
            echo "Error: 'generate' command requires <platform> <devices> <records>." >&2
            usage
        fi
        do_generate "$2" "$3" "$4"
        ;;

    run)
        if [ "$#" -ne 4 ]; then
            echo "Error: 'run' command requires <platform> <devices> <records>." >&2
            usage
        fi
        do_run "$2" "$3" "$4"
        ;;

    clean)
        echo "Cleaning up artifacts and data..."
        rm -f "$PROGRAM_NAME" "$MANIFEST" "$MANIFEST_SGX" "$DATA_FILE" "$KEY_FILE" "$TERSE_DATA_FILE" "$RESULTS_FILE"
        ;;

    *)
        echo "Error: Unknown command '$COMMAND'." >&2
        usage
        ;;
esac

echo "Done."
