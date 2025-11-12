#!/bin/bash

# aggregate_results.sh - Aggregate all benchmark results into a single file

OUTPUT_FILE="all_results.csv"

echo "Aggregating results into $OUTPUT_FILE..."

# Write header
echo "benchmark,platform,num_devices,records_per_device,avg_read_per_round_us,avg_sum_per_round_us,avg_lookup_per_round_us,avg_total_per_round_us" > "$OUTPUT_FILE"

# Process each benchmark directory
for dir in aes aes-shm terse terse-shm; do
    if [ -d "$dir" ]; then
        # Find the results CSV file
        results_file=$(find "$dir" -name "*_results.csv" -type f | head -n 1)

        if [ -f "$results_file" ]; then
            echo "Processing $results_file..."

            # Skip header and add benchmark name as first column
            tail -n +2 "$results_file" | while IFS=, read -r platform devices records read sum lookup total; do
                # For AES benchmarks, sum and lookup are combined in aggregation
                if [[ "$dir" == "aes"* ]]; then
                    # AES format: platform,devices,records,read,aggregation,write,total
                    # Map: aggregation->sum, write->lookup
                    echo "$dir,$platform,$devices,$records,$read,$sum,$lookup,$total"
                else
                    # Terse format already matches
                    echo "$dir,$platform,$devices,$records,$read,$sum,$lookup,$total"
                fi
            done >> "$OUTPUT_FILE"
        else
            echo "Warning: No results file found in $dir/"
        fi
    fi
done

echo "Aggregation complete: $OUTPUT_FILE"
echo ""
echo "Summary:"
wc -l "$OUTPUT_FILE"
echo ""
head -n 1 "$OUTPUT_FILE"
echo "..."
tail -n 5 "$OUTPUT_FILE"
