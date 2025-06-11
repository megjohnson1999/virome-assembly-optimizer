#!/bin/bash

# Individual assembly script for viral metagenomic samples
# Assembles each sample separately using both MEGAHIT and metaSPAdes

set -euo pipefail

# Configuration
INPUT_DIR=${INPUT_DIR:-"data/cleaned_reads"}
OUTPUT_DIR=${OUTPUT_DIR:-"results/assemblies/individual"}
THREADS=${THREADS:-8}
MEMORY=${MEMORY:-32}  # GB
MIN_CONTIG_LEN=${MIN_CONTIG_LEN:-500}
ASSEMBLERS=${ASSEMBLERS:-"megahit metaspades"}  # Space-separated list

# Create output directories
mkdir -p "$OUTPUT_DIR"/{megahit,metaspades,logs,stats}

echo "Starting individual assembly pipeline..."
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "Memory: ${MEMORY}GB"
echo "Assemblers: $ASSEMBLERS"
echo "Minimum contig length: $MIN_CONTIG_LEN"

# Function to run MEGAHIT assembly
run_megahit() {
    local sample_name=$1
    local r1_file=$2
    local r2_file=$3
    local output_dir=$4
    
    echo "  Running MEGAHIT for $sample_name..."
    
    local sample_output="$output_dir/megahit/${sample_name}"
    
    # Remove existing output directory if it exists
    if [[ -d "$sample_output" ]]; then
        rm -rf "$sample_output"
    fi
    
    # Run MEGAHIT
    megahit \
        -1 "$r1_file" \
        -2 "$r2_file" \
        -o "$sample_output" \
        --num-cpu-threads "$THREADS" \
        --memory 0.8 \
        --min-contig-len "$MIN_CONTIG_LEN" \
        --presets meta-sensitive \
        2>&1 | tee "$output_dir/logs/${sample_name}_megahit.log"
    
    # Rename final contigs
    if [[ -f "$sample_output/final.contigs.fa" ]]; then
        cp "$sample_output/final.contigs.fa" "$output_dir/megahit/${sample_name}_contigs.fasta"
        
        # Add sample prefix to contig names
        sed -i "s/^>/>${sample_name}_MEGAHIT_/" "$output_dir/megahit/${sample_name}_contigs.fasta"
        
        echo "    MEGAHIT assembly completed: $output_dir/megahit/${sample_name}_contigs.fasta"
    else
        echo "    ERROR: MEGAHIT assembly failed for $sample_name"
        return 1
    fi
}

# Function to run metaSPAdes assembly
run_metaspades() {
    local sample_name=$1
    local r1_file=$2
    local r2_file=$3
    local output_dir=$4
    
    echo "  Running metaSPAdes for $sample_name..."
    
    local sample_output="$output_dir/metaspades/${sample_name}"
    
    # Remove existing output directory if it exists
    if [[ -d "$sample_output" ]]; then
        rm -rf "$sample_output"
    fi
    
    # Run metaSPAdes
    metaspades.py \
        -1 "$r1_file" \
        -2 "$r2_file" \
        -o "$sample_output" \
        --threads "$THREADS" \
        --memory "$MEMORY" \
        2>&1 | tee "$output_dir/logs/${sample_name}_metaspades.log"
    
    # Copy and rename final contigs
    if [[ -f "$sample_output/contigs.fasta" ]]; then
        cp "$sample_output/contigs.fasta" "$output_dir/metaspades/${sample_name}_contigs.fasta"
        
        # Filter by minimum length and add sample prefix
        seqtk seq -L "$MIN_CONTIG_LEN" "$sample_output/contigs.fasta" | \
        sed "s/^>/>${sample_name}_SPADES_/" > "$output_dir/metaspades/${sample_name}_contigs.fasta"
        
        echo "    metaSPAdes assembly completed: $output_dir/metaspades/${sample_name}_contigs.fasta"
    else
        echo "    ERROR: metaSPAdes assembly failed for $sample_name"
        return 1
    fi
}

# Function to calculate assembly statistics
calculate_stats() {
    local sample_name=$1
    local assembler=$2
    local contigs_file=$3
    local output_dir=$4
    
    if [[ ! -f "$contigs_file" ]]; then
        echo "Warning: Contigs file not found: $contigs_file"
        return 1
    fi
    
    echo "  Calculating statistics for $sample_name ($assembler)..."
    
    # Use seqkit for statistics
    seqkit stats "$contigs_file" > "$output_dir/stats/${sample_name}_${assembler}_stats.txt"
    
    # Calculate N50 and other metrics using a simple awk script
    awk '
    BEGIN { 
        total_length = 0; 
        contig_count = 0; 
    }
    /^>/ { 
        if (length > 0) {
            lengths[++contig_count] = length;
            total_length += length;
        }
        length = 0; 
        next; 
    }
    { 
        length += length($0); 
    }
    END {
        if (length > 0) {
            lengths[++contig_count] = length;
            total_length += length;
        }
        
        # Sort lengths in descending order
        n = asort(lengths, sorted_lengths, "@val_num_desc");
        
        # Calculate N50
        target = total_length / 2;
        cumulative = 0;
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
        for (i = 1; i <= n; i++) {
            cumulative += sorted_lengths[i];
            if (cumulative >= target90) {
                n90 = sorted_lengths[i];
                break;
            }
        }
        
        printf "Sample: %s\n", "'$sample_name'";
        printf "Assembler: %s\n", "'$assembler'";
        printf "Contigs: %d\n", contig_count;
        printf "Total_length: %d\n", total_length;
        printf "Mean_length: %.2f\n", total_length / contig_count;
        printf "Max_length: %d\n", sorted_lengths[1];
        printf "Min_length: %d\n", sorted_lengths[n];
        printf "N50: %d\n", n50;
        printf "N90: %d\n", n90;
    }' "$contigs_file" >> "$output_dir/stats/${sample_name}_${assembler}_detailed_stats.txt"
}

# Main processing loop
total_samples=0
successful_assemblies=0
failed_assemblies=0

# Find all paired-end fastq files
find "$INPUT_DIR" -name "*_R1*.fastq*" | sort | while read -r r1_file; do
    # Extract sample name
    sample_name=$(basename "$r1_file" | sed 's/_R[12].*fastq.*//')
    r2_file=${r1_file/_R1/_R2}
    
    # Check if R2 file exists
    if [[ ! -f "$r2_file" ]]; then
        echo "Warning: R2 file not found for $sample_name, skipping..."
        continue
    fi
    
    echo "Processing sample: $sample_name"
    echo "  R1: $r1_file"
    echo "  R2: $r2_file"
    
    ((total_samples++))
    sample_success=true
    
    # Run assemblies based on configuration
    for assembler in $ASSEMBLERS; do
        case $assembler in
            "megahit")
                if ! run_megahit "$sample_name" "$r1_file" "$r2_file" "$OUTPUT_DIR"; then
                    sample_success=false
                else
                    calculate_stats "$sample_name" "megahit" "$OUTPUT_DIR/megahit/${sample_name}_contigs.fasta" "$OUTPUT_DIR"
                fi
                ;;
            "metaspades")
                if ! run_metaspades "$sample_name" "$r1_file" "$r2_file" "$OUTPUT_DIR"; then
                    sample_success=false
                else
                    calculate_stats "$sample_name" "metaspades" "$OUTPUT_DIR/metaspades/${sample_name}_contigs.fasta" "$OUTPUT_DIR"
                fi
                ;;
            *)
                echo "Warning: Unknown assembler '$assembler', skipping..."
                ;;
        esac
    done
    
    if $sample_success; then
        ((successful_assemblies++))
        echo "  Sample $sample_name completed successfully"
    else
        ((failed_assemblies++))
        echo "  Sample $sample_name had assembly failures"
    fi
    
    echo ""
done

# Generate summary report
echo "Individual assembly pipeline completed!"
echo "=================================="
echo "Total samples processed: $total_samples"
echo "Successful assemblies: $successful_assemblies"
echo "Failed assemblies: $failed_assemblies"

# Create combined statistics file
combined_stats_file="$OUTPUT_DIR/stats/all_samples_summary.txt"
echo "Creating combined statistics file: $combined_stats_file"

echo "Sample,Assembler,Contigs,Total_length,Mean_length,Max_length,Min_length,N50,N90" > "$combined_stats_file"

for stats_file in "$OUTPUT_DIR"/stats/*_detailed_stats.txt; do
    if [[ -f "$stats_file" ]]; then
        # Parse the detailed stats file and convert to CSV format
        awk '
        BEGIN { 
            sample = ""; assembler = ""; contigs = ""; total = ""; mean = ""; 
            max_len = ""; min_len = ""; n50 = ""; n90 = "";
        }
        /^Sample:/ { sample = $2 }
        /^Assembler:/ { assembler = $2 }
        /^Contigs:/ { contigs = $2 }
        /^Total_length:/ { total = $2 }
        /^Mean_length:/ { mean = $2 }
        /^Max_length:/ { max_len = $2 }
        /^Min_length:/ { min_len = $2 }
        /^N50:/ { n50 = $2 }
        /^N90:/ { n90 = $2 }
        END {
            if (sample != "") {
                printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n", 
                       sample, assembler, contigs, total, mean, max_len, min_len, n50, n90
            }
        }' "$stats_file" >> "$combined_stats_file"
    fi
done

echo "Combined statistics saved to: $combined_stats_file"

# Create a simple summary plot using R if available
if command -v Rscript &> /dev/null; then
    echo "Generating assembly summary plots..."
    
    cat > "$OUTPUT_DIR/stats/plot_assembly_stats.R" << 'EOF'
library(ggplot2)
library(dplyr)

# Read data
stats <- read.csv("all_samples_summary.txt")

# Create plots directory
dir.create("plots", showWarnings = FALSE)

# N50 comparison plot
p1 <- ggplot(stats, aes(x = Sample, y = N50, fill = Assembler)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "N50 Comparison Across Samples", y = "N50 (bp)")

ggsave("plots/n50_comparison.png", p1, width = 12, height = 6, dpi = 300)

# Contig count comparison
p2 <- ggplot(stats, aes(x = Sample, y = Contigs, fill = Assembler)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Contig Count Comparison", y = "Number of Contigs")

ggsave("plots/contig_count_comparison.png", p2, width = 12, height = 6, dpi = 300)

# Total assembly size comparison  
p3 <- ggplot(stats, aes(x = Sample, y = Total_length/1000000, fill = Assembler)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Total Assembly Size Comparison", y = "Total Length (Mbp)")

ggsave("plots/total_length_comparison.png", p3, width = 12, height = 6, dpi = 300)

cat("Plots saved to plots/ directory\n")
EOF

    cd "$OUTPUT_DIR/stats"
    Rscript plot_assembly_stats.R
    cd - > /dev/null
fi

echo ""
echo "Individual assembly results:"
echo "  Assembly files: $OUTPUT_DIR/{megahit,metaspades}/"
echo "  Statistics: $OUTPUT_DIR/stats/"
echo "  Logs: $OUTPUT_DIR/logs/"
echo ""
echo "Next steps:"
echo "  1. Run quality assessment: bash scripts/quality_assessment/01_run_checkv.sh"
echo "  2. Compare with co-assembly strategies"