#!/bin/bash

# Strategic co-assembly script based on similarity and variable analysis
# Executes co-assemblies defined by create_coassembly_groups.py

set -euo pipefail

# Configuration
GROUPS_DIR=${GROUPS_DIR:-"results/coassembly_groups"}
OUTPUT_DIR=${OUTPUT_DIR:-"results/assemblies/strategic_coassembly"}
THREADS=${THREADS:-8}
PARALLEL_JOBS=${PARALLEL_JOBS:-2}  # Number of assemblies to run in parallel

echo "Starting strategic co-assembly pipeline..."
echo "Groups directory: $GROUPS_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Threads per assembly: $THREADS"
echo "Parallel jobs: $PARALLEL_JOBS"

# Create output directories
mkdir -p "$OUTPUT_DIR"/{assemblies,logs,stats}

# Check if co-assembly groups exist
COMMANDS_DIR="$GROUPS_DIR/commands"
if [[ ! -d "$COMMANDS_DIR" ]]; then
    echo "Error: Commands directory not found: $COMMANDS_DIR"
    echo "Please run create_coassembly_groups.py first to generate assembly groups"
    exit 1
fi

# Count available command scripts
COMMAND_SCRIPTS=($(find "$COMMANDS_DIR" -name "group_*.sh" | sort))
TOTAL_GROUPS=${#COMMAND_SCRIPTS[@]}

if [[ $TOTAL_GROUPS -eq 0 ]]; then
    echo "Error: No assembly command scripts found in $COMMANDS_DIR"
    echo "Please run create_coassembly_groups.py first"
    exit 1
fi

echo "Found $TOTAL_GROUPS assembly groups to process"

# Function to run a single assembly command
run_assembly() {
    local script_path=$1
    local script_name=$(basename "$script_path")
    local group_id=${script_name%.sh}
    local log_file="$OUTPUT_DIR/logs/${group_id}.log"
    
    echo "Starting assembly: $group_id"
    echo "  Script: $script_path"
    echo "  Log: $log_file"
    
    # Create assembly-specific output directory
    local assembly_output="$OUTPUT_DIR/assemblies/$group_id"
    mkdir -p "$assembly_output"
    
    # Modify the script to use our output directory structure
    # Create a temporary script with corrected paths
    local temp_script="/tmp/${group_id}_modified.sh"
    sed "s|OUTPUT_DIR='results/assemblies/strategic_coassembly/${group_id}'|OUTPUT_DIR='$assembly_output'|g" \
        "$script_path" > "$temp_script"
    chmod +x "$temp_script"
    
    # Run the assembly
    local start_time=$(date +%s)
    
    if bash "$temp_script" > "$log_file" 2>&1; then
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        echo "  Completed successfully in ${duration}s"
        
        # Generate assembly statistics
        generate_assembly_stats "$group_id" "$assembly_output" "$log_file"
        
    else
        echo "  Failed - check log: $log_file"
        echo "FAILED" > "$OUTPUT_DIR/stats/${group_id}_status.txt"
    fi
    
    # Clean up temporary script
    rm -f "$temp_script"
}

# Function to generate assembly statistics
generate_assembly_stats() {
    local group_id=$1
    local assembly_dir=$2
    local log_file=$3
    
    local stats_file="$OUTPUT_DIR/stats/${group_id}_stats.txt"
    
    echo "Generating statistics for $group_id..."
    
    # Find the main contigs file
    local contigs_file=""
    if [[ -f "$assembly_dir/${group_id}_contigs.fasta" ]]; then
        contigs_file="$assembly_dir/${group_id}_contigs.fasta"
    elif [[ -f "$assembly_dir/megahit_coassembly/final.contigs.fa" ]]; then
        contigs_file="$assembly_dir/megahit_coassembly/final.contigs.fa"
    else
        # Look for any .fasta files
        contigs_file=$(find "$assembly_dir" -name "*.fasta" -o -name "*.fa" | head -1)
    fi
    
    if [[ -z "$contigs_file" || ! -f "$contigs_file" ]]; then
        echo "Warning: No contigs file found for $group_id"
        echo "STATUS: NO_CONTIGS" > "$stats_file"
        return 1
    fi
    
    # Calculate basic statistics using awk
    awk '
    BEGIN { 
        total_length = 0; 
        contig_count = 0; 
        current_length = 0;
    }
    /^>/ { 
        if (current_length > 0) {
            lengths[++contig_count] = current_length;
            total_length += current_length;
        }
        current_length = 0; 
        next; 
    }
    { 
        current_length += length($0); 
    }
    END {
        if (current_length > 0) {
            lengths[++contig_count] = current_length;
            total_length += current_length;
        }
        
        if (contig_count == 0) {
            print "STATUS: EMPTY_FILE";
            exit;
        }
        
        # Sort lengths in descending order
        n = asort(lengths, sorted_lengths, "@val_num_desc");
        
        # Calculate N50
        target = total_length / 2;
        cumulative = 0;
        n50 = 0;
        for (i = 1; i <= n; i++) {
            cumulative += sorted_lengths[i];
            if (cumulative >= target) {
                n50 = sorted_lengths[i];
                break;
            }
        }
        
        # Calculate N90
        target90 = total_length * 0.9;
        cumulative = 0;
        n90 = 0;
        for (i = 1; i <= n; i++) {
            cumulative += sorted_lengths[i];
            if (cumulative >= target90) {
                n90 = sorted_lengths[i];
                break;
            }
        }
        
        printf "GROUP_ID: %s\n", "'$group_id'";
        printf "STATUS: SUCCESS\n";
        printf "CONTIGS_FILE: %s\n", "'$contigs_file'";
        printf "CONTIGS_COUNT: %d\n", contig_count;
        printf "TOTAL_LENGTH: %d\n", total_length;
        printf "MEAN_LENGTH: %.2f\n", total_length / contig_count;
        printf "MAX_LENGTH: %d\n", sorted_lengths[1];
        printf "MIN_LENGTH: %d\n", sorted_lengths[n];
        printf "N50: %d\n", n50;
        printf "N90: %d\n", n90;
        
        # Categorize contigs by length
        short_contigs = 0;  # < 1kb
        medium_contigs = 0; # 1kb - 10kb  
        long_contigs = 0;   # > 10kb
        
        for (i = 1; i <= n; i++) {
            if (sorted_lengths[i] < 1000) {
                short_contigs++;
            } else if (sorted_lengths[i] < 10000) {
                medium_contigs++;
            } else {
                long_contigs++;
            }
        }
        
        printf "SHORT_CONTIGS: %d\n", short_contigs;
        printf "MEDIUM_CONTIGS: %d\n", medium_contigs;
        printf "LONG_CONTIGS: %d\n", long_contigs;
    }' "$contigs_file" > "$stats_file"
    
    # Add timing information from log
    if [[ -f "$log_file" ]]; then
        echo "" >> "$stats_file"
        echo "TIMING_INFO:" >> "$stats_file"
        
        # Extract MEGAHIT timing if available
        if grep -q "MEGAHIT" "$log_file"; then
            echo "ASSEMBLER: MEGAHIT" >> "$stats_file"
            
            # Try to extract timing information
            local start_time=$(grep "Start MEGAHIT" "$log_file" | head -1 | cut -d' ' -f1-2 2>/dev/null || echo "")
            local end_time=$(grep "ALL DONE" "$log_file" | head -1 | cut -d' ' -f1-2 2>/dev/null || echo "")
            
            if [[ -n "$start_time" && -n "$end_time" ]]; then
                echo "START_TIME: $start_time" >> "$stats_file"
                echo "END_TIME: $end_time" >> "$stats_file"
            fi
        fi
        
        # Add memory usage if available
        local max_memory=$(grep -o "maximum resident set size.*" "$log_file" | tail -1 || echo "")
        if [[ -n "$max_memory" ]]; then
            echo "MAX_MEMORY: $max_memory" >> "$stats_file"
        fi
    fi
}

# Function to run assemblies in parallel
run_assemblies_parallel() {
    local scripts=("$@")
    local total=${#scripts[@]}
    local completed=0
    local failed=0
    
    echo "Running $total assemblies with up to $PARALLEL_JOBS parallel jobs..."
    
    # Use GNU parallel if available, otherwise use xargs
    if command -v parallel &> /dev/null; then
        echo "Using GNU parallel for job execution"
        printf '%s\n' "${scripts[@]}" | \
        parallel -j "$PARALLEL_JOBS" --bar run_assembly
    else
        echo "Using xargs for job execution"
        printf '%s\n' "${scripts[@]}" | \
        xargs -I {} -P "$PARALLEL_JOBS" bash -c 'run_assembly "$@"' _ {}
    fi
}

# Export the function so it can be used by parallel/xargs
export -f run_assembly
export -f generate_assembly_stats
export OUTPUT_DIR
export THREADS

# Main execution
echo "Starting strategic co-assembly execution..."

# Option 1: Run assemblies in parallel
if [[ $PARALLEL_JOBS -gt 1 ]]; then
    run_assemblies_parallel "${COMMAND_SCRIPTS[@]}"
else
    # Option 2: Run assemblies sequentially
    for script in "${COMMAND_SCRIPTS[@]}"; do
        run_assembly "$script"
    done
fi

# Generate final summary
echo "Generating final summary..."

SUMMARY_FILE="$OUTPUT_DIR/strategic_coassembly_summary.txt"
CSV_FILE="$OUTPUT_DIR/strategic_coassembly_summary.csv"

# Create text summary
cat > "$SUMMARY_FILE" << EOF
Strategic Co-assembly Pipeline Summary
======================================

Execution Date: $(date)
Total Groups: $TOTAL_GROUPS
Parallel Jobs: $PARALLEL_JOBS
Threads per Job: $THREADS

Group Results:
EOF

# Create CSV header
echo "Group_ID,Status,Contigs_Count,Total_Length,Mean_Length,Max_Length,N50,N90,Short_Contigs,Medium_Contigs,Long_Contigs" > "$CSV_FILE"

# Process results
successful_groups=0
failed_groups=0

for script in "${COMMAND_SCRIPTS[@]}"; do
    script_name=$(basename "$script")
    group_id=${script_name%.sh}
    stats_file="$OUTPUT_DIR/stats/${group_id}_stats.txt"
    
    echo "" >> "$SUMMARY_FILE"
    echo "Group: $group_id" >> "$SUMMARY_FILE"
    
    if [[ -f "$stats_file" ]]; then
        # Read status
        status=$(grep "^STATUS:" "$stats_file" | cut -d' ' -f2 || echo "UNKNOWN")
        
        if [[ "$status" == "SUCCESS" ]]; then
            ((successful_groups++))
            echo "  Status: SUCCESS" >> "$SUMMARY_FILE"
            
            # Extract key metrics for CSV
            contigs=$(grep "^CONTIGS_COUNT:" "$stats_file" | cut -d' ' -f2 || echo "0")
            total_len=$(grep "^TOTAL_LENGTH:" "$stats_file" | cut -d' ' -f2 || echo "0")
            mean_len=$(grep "^MEAN_LENGTH:" "$stats_file" | cut -d' ' -f2 || echo "0")
            max_len=$(grep "^MAX_LENGTH:" "$stats_file" | cut -d' ' -f2 || echo "0")
            n50=$(grep "^N50:" "$stats_file" | cut -d' ' -f2 || echo "0")
            n90=$(grep "^N90:" "$stats_file" | cut -d' ' -f2 || echo "0")
            short_contigs=$(grep "^SHORT_CONTIGS:" "$stats_file" | cut -d' ' -f2 || echo "0")
            medium_contigs=$(grep "^MEDIUM_CONTIGS:" "$stats_file" | cut -d' ' -f2 || echo "0")
            long_contigs=$(grep "^LONG_CONTIGS:" "$stats_file" | cut -d' ' -f2 || echo "0")
            
            echo "  Contigs: $contigs" >> "$SUMMARY_FILE"
            echo "  Total Length: $total_len bp" >> "$SUMMARY_FILE"
            echo "  N50: $n50 bp" >> "$SUMMARY_FILE"
            
            # Add to CSV
            echo "$group_id,$status,$contigs,$total_len,$mean_len,$max_len,$n50,$n90,$short_contigs,$medium_contigs,$long_contigs" >> "$CSV_FILE"
        else
            ((failed_groups++))
            echo "  Status: $status" >> "$SUMMARY_FILE"
            echo "$group_id,$status,0,0,0,0,0,0,0,0,0" >> "$CSV_FILE"
        fi
    else
        ((failed_groups++))
        echo "  Status: NO_STATS_FILE" >> "$SUMMARY_FILE"
        echo "$group_id,NO_STATS_FILE,0,0,0,0,0,0,0,0,0" >> "$CSV_FILE"
    fi
done

# Add final statistics to summary
cat >> "$SUMMARY_FILE" << EOF

Final Statistics:
=================
Successful Groups: $successful_groups
Failed Groups: $failed_groups
Success Rate: $(echo "scale=1; $successful_groups * 100 / $TOTAL_GROUPS" | bc -l)%

Output Files:
- Assemblies: $OUTPUT_DIR/assemblies/
- Statistics: $OUTPUT_DIR/stats/
- Logs: $OUTPUT_DIR/logs/
- Summary: $SUMMARY_FILE
- CSV Data: $CSV_FILE
EOF

echo ""
echo "Strategic co-assembly pipeline completed!"
echo "=================================="
echo "Total groups: $TOTAL_GROUPS"
echo "Successful: $successful_groups"
echo "Failed: $failed_groups"
echo "Success rate: $(echo "scale=1; $successful_groups * 100 / $TOTAL_GROUPS" | bc -l)%"
echo ""
echo "Results saved to:"
echo "  Summary: $SUMMARY_FILE"
echo "  CSV data: $CSV_FILE"
echo "  Assemblies: $OUTPUT_DIR/assemblies/"
echo ""
echo "Next steps:"
echo "  1. Run quality assessment: bash scripts/quality_assessment/01_run_checkv.sh"
echo "  2. Compare with other assembly strategies"