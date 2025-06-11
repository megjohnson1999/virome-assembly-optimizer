#!/bin/bash

# Global co-assembly script - assembles all samples together
# Useful as a comparison strategy and for highly similar sample sets

set -euo pipefail

# Configuration
INPUT_DIR=${INPUT_DIR:-"data/cleaned_reads"}
OUTPUT_DIR=${OUTPUT_DIR:-"results/assemblies/global_coassembly"}
THREADS=${THREADS:-16}
MEMORY=${MEMORY:-64}  # GB - Global assemblies need more memory
MIN_CONTIG_LEN=${MIN_CONTIG_LEN:-500}
ASSEMBLER=${ASSEMBLER:-"megahit"}  # megahit or metaspades
MAX_SAMPLES=${MAX_SAMPLES:-50}  # Safety limit for computational resources

echo "Starting global co-assembly pipeline..."
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Assembler: $ASSEMBLER"
echo "Threads: $THREADS"
echo "Memory: ${MEMORY}GB"
echo "Minimum contig length: $MIN_CONTIG_LEN"
echo "Maximum samples: $MAX_SAMPLES"

# Create output directories
mkdir -p "$OUTPUT_DIR"/{assembly,logs,stats,temp}

# Function to find and validate input files
find_input_files() {
    echo "Discovering input files..."
    
    local r1_files=()
    local r2_files=()
    local sample_names=()
    
    # Find all R1 files
    while IFS= read -r -d '' r1_file; do
        # Extract sample name
        local sample_name=$(basename "$r1_file" | sed 's/_R[12].*fastq.*//')
        local r2_file=${r1_file/_R1/_R2}
        
        # Check if R2 file exists
        if [[ -f "$r2_file" ]]; then
            r1_files+=("$r1_file")
            r2_files+=("$r2_file")
            sample_names+=("$sample_name")
            echo "  Found: $sample_name"
        else
            echo "  Warning: R2 file missing for $sample_name, skipping"
        fi
    done < <(find "$INPUT_DIR" -name "*_R1*.fastq*" -print0 | sort -z)
    
    local n_samples=${#sample_names[@]}
    echo "Found $n_samples valid sample pairs"
    
    if [[ $n_samples -eq 0 ]]; then
        echo "Error: No valid paired-end samples found in $INPUT_DIR"
        exit 1
    fi
    
    if [[ $n_samples -gt $MAX_SAMPLES ]]; then
        echo "Warning: Found $n_samples samples, but limit is $MAX_SAMPLES"
        echo "Consider using strategic co-assembly instead for better resource management"
        echo "Proceeding with first $MAX_SAMPLES samples..."
        
        # Truncate arrays
        r1_files=("${r1_files[@]:0:$MAX_SAMPLES}")
        r2_files=("${r2_files[@]:0:$MAX_SAMPLES}")
        sample_names=("${sample_names[@]:0:$MAX_SAMPLES}")
        n_samples=$MAX_SAMPLES
    fi
    
    # Estimate resource requirements
    local avg_file_size=$(du -sb "${r1_files[@]}" | awk '{sum+=$1} END {print int(sum/NR)}')
    local total_data_gb=$((avg_file_size * n_samples * 2 / 1024 / 1024 / 1024))
    local estimated_memory_gb=$((total_data_gb / 4))  # Rough estimate
    
    echo ""
    echo "Resource Estimation:"
    echo "  Samples: $n_samples"
    echo "  Total data: ~${total_data_gb}GB"
    echo "  Estimated memory needed: ~${estimated_memory_gb}GB"
    echo "  Configured memory: ${MEMORY}GB"
    
    if [[ $estimated_memory_gb -gt $MEMORY ]]; then
        echo "  Warning: May need more memory than configured"
    fi
    
    # Save file lists for assembly
    printf '%s\n' "${r1_files[@]}" > "$OUTPUT_DIR/temp/r1_files.txt"
    printf '%s\n' "${r2_files[@]}" > "$OUTPUT_DIR/temp/r2_files.txt"
    printf '%s\n' "${sample_names[@]}" > "$OUTPUT_DIR/temp/sample_names.txt"
    
    echo "Input files saved to: $OUTPUT_DIR/temp/"
    echo ""
}

# Function to run MEGAHIT global assembly
run_megahit_global() {
    echo "Running MEGAHIT global co-assembly..."
    
    local r1_list=$(cat "$OUTPUT_DIR/temp/r1_files.txt" | tr '\n' ',')
    local r2_list=$(cat "$OUTPUT_DIR/temp/r2_files.txt" | tr '\n' ',')
    
    # Remove trailing comma
    r1_list=${r1_list%,}
    r2_list=${r2_list%,}
    
    local assembly_dir="$OUTPUT_DIR/assembly/megahit_global"
    local log_file="$OUTPUT_DIR/logs/megahit_global.log"
    
    echo "  Output directory: $assembly_dir"
    echo "  Log file: $log_file"
    
    # Remove existing output directory if it exists
    if [[ -d "$assembly_dir" ]]; then
        echo "  Removing existing assembly directory..."
        rm -rf "$assembly_dir"
    fi
    
    # Run MEGAHIT
    echo "  Starting MEGAHIT assembly..."
    local start_time=$(date +%s)
    
    megahit \
        -1 "$r1_list" \
        -2 "$r2_list" \
        -o "$assembly_dir" \
        --num-cpu-threads "$THREADS" \
        --memory 0.8 \
        --min-contig-len "$MIN_CONTIG_LEN" \
        --presets meta-sensitive \
        --verbose \
        2>&1 | tee "$log_file"
    
    local exit_code=${PIPESTATUS[0]}
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    if [[ $exit_code -eq 0 ]]; then
        echo "  MEGAHIT completed successfully in ${duration}s"
        
        # Copy and rename final contigs
        if [[ -f "$assembly_dir/final.contigs.fa" ]]; then
            cp "$assembly_dir/final.contigs.fa" "$OUTPUT_DIR/global_coassembly_contigs.fasta"
            
            # Add prefix to contig names
            sed -i 's/^>/>\GLOBAL_COASM_/' "$OUTPUT_DIR/global_coassembly_contigs.fasta"
            
            echo "  Final contigs: $OUTPUT_DIR/global_coassembly_contigs.fasta"
            return 0
        else
            echo "  Error: Final contigs file not found"
            return 1
        fi
    else
        echo "  Error: MEGAHIT failed with exit code $exit_code"
        return 1
    fi
}

# Function to run metaSPAdes global assembly
run_metaspades_global() {
    echo "Running metaSPAdes global co-assembly..."
    
    local assembly_dir="$OUTPUT_DIR/assembly/metaspades_global"
    local log_file="$OUTPUT_DIR/logs/metaspades_global.log"
    
    echo "  Output directory: $assembly_dir"
    echo "  Log file: $log_file"
    
    # Remove existing output directory if it exists
    if [[ -d "$assembly_dir" ]]; then
        echo "  Removing existing assembly directory..."
        rm -rf "$assembly_dir"
    fi
    
    # Prepare file list arguments
    local r1_args=""
    local r2_args=""
    
    while IFS= read -r r1_file; do
        r1_args="$r1_args $r1_file"
    done < "$OUTPUT_DIR/temp/r1_files.txt"
    
    while IFS= read -r r2_file; do
        r2_args="$r2_args $r2_file"
    done < "$OUTPUT_DIR/temp/r2_files.txt"
    
    # Run metaSPAdes
    echo "  Starting metaSPAdes assembly..."
    local start_time=$(date +%s)
    
    metaspades.py \
        -1 $r1_args \
        -2 $r2_args \
        -o "$assembly_dir" \
        --threads "$THREADS" \
        --memory "$MEMORY" \
        2>&1 | tee "$log_file"
    
    local exit_code=${PIPESTATUS[0]}
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    if [[ $exit_code -eq 0 ]]; then
        echo "  metaSPAdes completed successfully in ${duration}s"
        
        # Filter and rename contigs
        if [[ -f "$assembly_dir/contigs.fasta" ]]; then
            seqtk seq -L "$MIN_CONTIG_LEN" "$assembly_dir/contigs.fasta" | \
            sed 's/^>/>\GLOBAL_SPADES_/' > "$OUTPUT_DIR/global_coassembly_contigs.fasta"
            
            echo "  Final contigs: $OUTPUT_DIR/global_coassembly_contigs.fasta"
            return 0
        else
            echo "  Error: Contigs file not found"
            return 1
        fi
    else
        echo "  Error: metaSPAdes failed with exit code $exit_code"
        return 1
    fi
}

# Function to calculate assembly statistics
calculate_assembly_stats() {
    local contigs_file="$OUTPUT_DIR/global_coassembly_contigs.fasta"
    local stats_file="$OUTPUT_DIR/stats/global_assembly_stats.txt"
    
    if [[ ! -f "$contigs_file" ]]; then
        echo "Error: Contigs file not found: $contigs_file"
        return 1
    fi
    
    echo "Calculating assembly statistics..."
    
    # Calculate comprehensive statistics using awk
    awk '
    BEGIN { 
        total_length = 0; 
        contig_count = 0; 
        current_length = 0;
        gc_content = 0;
        gc_count = 0;
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
        seq = toupper($0);
        current_length += length(seq);
        
        # Count GC content
        for (i = 1; i <= length(seq); i++) {
            char = substr(seq, i, 1);
            if (char == "G" || char == "C") {
                gc_count++;
            }
        }
    }
    END {
        if (current_length > 0) {
            lengths[++contig_count] = current_length;
            total_length += current_length;
        }
        
        if (contig_count == 0) {
            print "ERROR: No contigs found";
            exit 1;
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
        
        # Calculate L50 (number of contigs needed to reach N50)
        cumulative = 0;
        l50 = 0;
        for (i = 1; i <= n; i++) {
            cumulative += sorted_lengths[i];
            l50 = i;
            if (cumulative >= target) {
                break;
            }
        }
        
        # Calculate GC content
        gc_percent = (total_length > 0) ? (gc_count / total_length) * 100 : 0;
        
        # Categorize contigs by length
        short_contigs = 0;   # < 1kb
        medium_contigs = 0;  # 1kb - 10kb  
        long_contigs = 0;    # 10kb - 100kb
        very_long_contigs = 0; # > 100kb
        
        for (i = 1; i <= n; i++) {
            if (sorted_lengths[i] < 1000) {
                short_contigs++;
            } else if (sorted_lengths[i] < 10000) {
                medium_contigs++;
            } else if (sorted_lengths[i] < 100000) {
                long_contigs++;
            } else {
                very_long_contigs++;
            }
        }
        
        # Print results
        print "GLOBAL_COASSEMBLY_STATISTICS";
        print "============================";
        print "";
        print "Basic Statistics:";
        print "  Assembler: " "'$ASSEMBLER'";
        print "  Total Contigs: " contig_count;
        print "  Total Length: " total_length " bp";
        print "  Mean Length: " sprintf("%.2f", total_length / contig_count) " bp";
        print "  Median Length: " sorted_lengths[int(n/2)] " bp";
        print "  Max Length: " sorted_lengths[1] " bp";
        print "  Min Length: " sorted_lengths[n] " bp";
        print "  GC Content: " sprintf("%.2f", gc_percent) "%";
        print "";
        print "Assembly Metrics:";
        print "  N50: " n50 " bp";
        print "  L50: " l50 " contigs";
        print "  N90: " n90 " bp";
        print "";
        print "Contig Length Distribution:";
        print "  Very Long (>100kb): " very_long_contigs " (" sprintf("%.1f", very_long_contigs/contig_count*100) "%)";
        print "  Long (10-100kb): " long_contigs " (" sprintf("%.1f", long_contigs/contig_count*100) "%)";
        print "  Medium (1-10kb): " medium_contigs " (" sprintf("%.1f", medium_contigs/contig_count*100) "%)";
        print "  Short (<1kb): " short_contigs " (" sprintf("%.1f", short_contigs/contig_count*100) "%)";
        print "";
        
        # CSV format for easy parsing
        print "CSV_FORMAT:";
        print "Metric,Value";
        print "Assembler," "'$ASSEMBLER'";
        print "Total_Contigs," contig_count;
        print "Total_Length," total_length;
        print "Mean_Length," sprintf("%.2f", total_length / contig_count);
        print "Max_Length," sorted_lengths[1];
        print "Min_Length," sorted_lengths[n];
        print "N50," n50;
        print "L50," l50;
        print "N90," n90;
        print "GC_Content," sprintf("%.2f", gc_percent);
        print "Very_Long_Contigs," very_long_contigs;
        print "Long_Contigs," long_contigs;
        print "Medium_Contigs," medium_contigs;
        print "Short_Contigs," short_contigs;
    }' "$contigs_file" > "$stats_file"
    
    echo "Statistics saved to: $stats_file"
    
    # Also create a simple CSV file
    local csv_file="$OUTPUT_DIR/stats/global_assembly_summary.csv"
    awk '/^CSV_FORMAT:/,0' "$stats_file" | tail -n +2 > "$csv_file"
    echo "CSV summary saved to: $csv_file"
    
    # Display key statistics
    echo ""
    echo "Key Assembly Statistics:"
    echo "======================="
    grep -E "(Total Contigs|Total Length|N50|GC Content)" "$stats_file" | sed 's/^  //'
}

# Function to create comparison with sample count
create_sample_info() {
    local info_file="$OUTPUT_DIR/sample_info.txt"
    local n_samples=$(wc -l < "$OUTPUT_DIR/temp/sample_names.txt")
    
    cat > "$info_file" << EOF
Global Co-assembly Information
==============================

Assembly Date: $(date)
Assembler: $ASSEMBLER
Input Directory: $INPUT_DIR
Number of Samples: $n_samples
Threads Used: $THREADS
Memory Allocated: ${MEMORY}GB

Samples Included:
EOF
    
    cat "$OUTPUT_DIR/temp/sample_names.txt" | nl -w3 -s'. ' >> "$info_file"
    
    echo ""
    echo "Sample information saved to: $info_file"
}

# Main execution
echo "Starting global co-assembly workflow..."
echo ""

# Step 1: Find and validate input files
find_input_files

# Step 2: Run assembly based on selected assembler
case $ASSEMBLER in
    "megahit")
        if run_megahit_global; then
            assembly_success=true
        else
            assembly_success=false
        fi
        ;;
    "metaspades")
        if run_metaspades_global; then
            assembly_success=true
        else
            assembly_success=false
        fi
        ;;
    *)
        echo "Error: Unknown assembler '$ASSEMBLER'"
        echo "Supported assemblers: megahit, metaspades"
        exit 1
        ;;
esac

# Step 3: Process results if assembly was successful
if [[ $assembly_success == true ]]; then
    echo ""
    echo "Assembly completed successfully!"
    
    # Calculate statistics
    calculate_assembly_stats
    
    # Create sample information
    create_sample_info
    
    # Clean up temporary files
    echo ""
    echo "Cleaning up temporary files..."
    rm -rf "$OUTPUT_DIR/temp"
    
    echo ""
    echo "Global co-assembly pipeline completed successfully!"
    echo "=================================================="
    echo ""
    echo "Output files:"
    echo "  Contigs: $OUTPUT_DIR/global_coassembly_contigs.fasta"
    echo "  Statistics: $OUTPUT_DIR/stats/global_assembly_stats.txt"
    echo "  CSV Summary: $OUTPUT_DIR/stats/global_assembly_summary.csv"
    echo "  Sample Info: $OUTPUT_DIR/sample_info.txt"
    echo "  Logs: $OUTPUT_DIR/logs/"
    echo ""
    echo "Next steps:"
    echo "  1. Run quality assessment: bash scripts/quality_assessment/01_run_checkv.sh"
    echo "  2. Compare with individual and strategic co-assembly results"
    echo "  3. Perform read mapping analysis"
    
else
    echo ""
    echo "Global co-assembly failed!"
    echo "========================="
    echo ""
    echo "Check the log files in $OUTPUT_DIR/logs/ for error details"
    echo ""
    echo "Common solutions:"
    echo "  1. Increase memory allocation (current: ${MEMORY}GB)"
    echo "  2. Reduce number of samples (current: $(wc -l < "$OUTPUT_DIR/temp/sample_names.txt" 2>/dev/null || echo 'unknown'))"
    echo "  3. Check input file quality and format"
    echo "  4. Try using MEGAHIT instead of metaSPAdes for lower memory usage"
    
    exit 1
fi