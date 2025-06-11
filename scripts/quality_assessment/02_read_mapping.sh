#!/bin/bash

# Read mapping analysis for assembly quality assessment
# Maps reads back to contigs to evaluate coverage, mapping rates, and assembly accuracy

set -euo pipefail

# Configuration
ASSEMBLIES_DIR=${ASSEMBLIES_DIR:-"results/assemblies"}
READS_DIR=${READS_DIR:-"data/cleaned_reads"}
OUTPUT_DIR=${OUTPUT_DIR:-"results/quality_assessment/read_mapping"}
THREADS=${THREADS:-8}
MIN_MAPQ=${MIN_MAPQ:-10}  # Minimum mapping quality
MIN_CONTIG_LEN=${MIN_CONTIG_LEN:-1000}  # Only analyze longer contigs

echo "Starting read mapping quality assessment..."
echo "Assemblies directory: $ASSEMBLIES_DIR"
echo "Reads directory: $READS_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "Minimum mapping quality: $MIN_MAPQ"
echo "Minimum contig length: $MIN_CONTIG_LEN"

# Create output directories
mkdir -p "$OUTPUT_DIR"/{individual,strategic_coassembly,global_coassembly,meta_assembly,combined,logs,temp}

# Function to check dependencies
check_dependencies() {
    echo "Checking dependencies..."
    
    local missing_deps=()
    
    if ! command -v bwa &> /dev/null; then
        missing_deps+=("bwa")
    fi
    
    if ! command -v samtools &> /dev/null; then
        missing_deps+=("samtools")
    fi
    
    if ! command -v bedtools &> /dev/null; then
        missing_deps+=("bedtools")
    fi
    
    if [[ ${#missing_deps[@]} -gt 0 ]]; then
        echo "Error: Missing required dependencies:"
        printf '  %s\n' "${missing_deps[@]}"
        echo ""
        echo "Installation suggestions:"
        echo "  conda install -c bioconda bwa samtools bedtools"
        exit 1
    fi
    
    echo "All dependencies found"
}

# Function to find and validate read files
find_read_files() {
    echo "Discovering read files..."
    
    local read_pairs=()
    
    # Find all R1 files and check for corresponding R2 files
    while IFS= read -r -d '' r1_file; do
        local sample_name=$(basename "$r1_file" | sed 's/_R[12].*fastq.*//')
        local r2_file=${r1_file/_R1/_R2}
        
        if [[ -f "$r2_file" ]]; then
            read_pairs+=("$sample_name:$r1_file:$r2_file")
            echo "  Found: $sample_name"
        else
            echo "  Warning: R2 file missing for $sample_name"
        fi
    done < <(find "$READS_DIR" -name "*_R1*.fastq*" -print0 | sort -z)
    
    local n_samples=${#read_pairs[@]}
    echo "Found $n_samples valid read pairs"
    
    if [[ $n_samples -eq 0 ]]; then
        echo "Error: No valid paired-end read files found in $READS_DIR"
        exit 1
    fi
    
    # Save read pairs list
    printf '%s\n' "${read_pairs[@]}" > "$OUTPUT_DIR/temp/read_pairs.txt"
    
    echo "Read pairs saved to: $OUTPUT_DIR/temp/read_pairs.txt"
}

# Function to map reads to contigs and calculate statistics
map_reads_to_assembly() {
    local contigs_file=$1
    local strategy_name=$2
    local output_base=$3
    
    echo "Mapping reads to $strategy_name assembly..."
    echo "  Contigs: $contigs_file"
    echo "  Output: $output_base"
    
    if [[ ! -f "$contigs_file" ]]; then
        echo "  Warning: Contigs file not found: $contigs_file"
        return 1
    fi
    
    # Filter contigs by minimum length
    local filtered_contigs="$output_base/filtered_contigs.fasta"
    seqtk seq -L "$MIN_CONTIG_LEN" "$contigs_file" > "$filtered_contigs"
    
    local n_contigs=$(grep -c "^>" "$filtered_contigs" 2>/dev/null || echo "0")
    
    if [[ $n_contigs -eq 0 ]]; then
        echo "  Warning: No contigs found after length filtering"
        return 1
    fi
    
    echo "  Processing $n_contigs contigs ≥ $MIN_CONTIG_LEN bp"
    
    # Build BWA index
    echo "  Building BWA index..."
    local index_base="$output_base/contigs_index"
    
    bwa index -p "$index_base" "$filtered_contigs" 2>"$OUTPUT_DIR/logs/${strategy_name}_index.log"
    
    # Initialize summary statistics
    local stats_file="$output_base/mapping_statistics.txt"
    local csv_file="$output_base/mapping_summary.csv"
    
    cat > "$stats_file" << EOF
Read Mapping Statistics - $strategy_name
========================================

Assembly: $contigs_file
Strategy: $strategy_name
Analysis Date: $(date)
Contigs Analyzed: $n_contigs (≥ $MIN_CONTIG_LEN bp)

Per-Sample Mapping Results:
EOF
    
    echo "Sample,Total_Reads,Mapped_Reads,Mapping_Rate,Properly_Paired,Mean_Coverage,Median_Coverage,Covered_Bases,Coverage_Breadth" > "$csv_file"
    
    # Map reads from each sample
    local total_samples=0
    local total_reads=0
    local total_mapped=0
    
    while IFS=: read -r sample_name r1_file r2_file; do
        echo "    Mapping sample: $sample_name"
        
        local sample_bam="$output_base/${sample_name}.bam"
        local sample_stats="$output_base/${sample_name}_stats.txt"
        
        # Map reads with BWA MEM
        bwa mem -t "$THREADS" "$index_base" "$r1_file" "$r2_file" 2>>"$OUTPUT_DIR/logs/${strategy_name}_mapping.log" | \
        samtools view -@ "$THREADS" -q "$MIN_MAPQ" -f 2 -Sb - | \
        samtools sort -@ "$THREADS" -o "$sample_bam" -
        
        # Index BAM file
        samtools index "$sample_bam"
        
        # Calculate mapping statistics
        local read_count=$(samtools view -@ "$THREADS" -c "$sample_bam")
        local flagstat_output=$(samtools flagstat "$sample_bam")
        
        # Parse flagstat output
        local total_reads_sample=$(echo "$flagstat_output" | grep "in total" | cut -d' ' -f1)
        local mapped_reads=$(echo "$flagstat_output" | grep "mapped (" | head -1 | cut -d' ' -f1)
        local properly_paired=$(echo "$flagstat_output" | grep "properly paired" | cut -d' ' -f1)
        
        # Calculate mapping rate
        local mapping_rate=0
        if [[ $total_reads_sample -gt 0 ]]; then
            mapping_rate=$(echo "scale=2; $mapped_reads * 100 / $total_reads_sample" | bc -l)
        fi
        
        # Calculate coverage statistics
        local coverage_stats=$(calculate_coverage_stats "$sample_bam" "$filtered_contigs" "$sample_name")
        
        # Extract coverage metrics
        local mean_coverage=$(echo "$coverage_stats" | grep "Mean coverage:" | cut -d':' -f2 | xargs)
        local median_coverage=$(echo "$coverage_stats" | grep "Median coverage:" | cut -d':' -f2 | xargs)
        local covered_bases=$(echo "$coverage_stats" | grep "Covered bases:" | cut -d':' -f2 | xargs)
        local coverage_breadth=$(echo "$coverage_stats" | grep "Coverage breadth:" | cut -d':' -f2 | xargs)
        
        # Update summary statistics
        cat >> "$stats_file" << EOF

Sample: $sample_name
  Total reads: $total_reads_sample
  Mapped reads: $mapped_reads
  Mapping rate: $mapping_rate%
  Properly paired: $properly_paired
  Mean coverage: $mean_coverage
  Median coverage: $median_coverage
  Coverage breadth: $coverage_breadth%
EOF
        
        # Add to CSV
        echo "$sample_name,$total_reads_sample,$mapped_reads,$mapping_rate,$properly_paired,$mean_coverage,$median_coverage,$covered_bases,$coverage_breadth" >> "$csv_file"
        
        # Update totals
        ((total_samples++))
        total_reads=$((total_reads + total_reads_sample))
        total_mapped=$((total_mapped + mapped_reads))
        
    done < "$OUTPUT_DIR/temp/read_pairs.txt"
    
    # Calculate overall statistics
    local overall_mapping_rate=0
    if [[ $total_reads -gt 0 ]]; then
        overall_mapping_rate=$(echo "scale=2; $total_mapped * 100 / $total_reads" | bc -l)
    fi
    
    cat >> "$stats_file" << EOF

Overall Summary:
================
Total samples: $total_samples
Total reads: $total_reads
Total mapped reads: $total_mapped
Overall mapping rate: $overall_mapping_rate%

Files Generated:
  BAM files: $output_base/*.bam
  Statistics: $stats_file
  CSV summary: $csv_file
EOF
    
    echo "  Mapping completed for $strategy_name"
    echo "    Overall mapping rate: $overall_mapping_rate%"
    echo "    Results saved to: $output_base"
}

# Function to calculate coverage statistics
calculate_coverage_stats() {
    local bam_file=$1
    local contigs_file=$2
    local sample_name=$3
    
    # Calculate per-base coverage
    local coverage_file="${bam_file%.bam}_coverage.txt"
    
    bedtools genomecov -ibam "$bam_file" -d > "$coverage_file"
    
    # Calculate statistics using awk
    awk '
    BEGIN {
        total_bases = 0;
        covered_bases = 0;
        sum_coverage = 0;
        split("", coverages);
    }
    {
        total_bases++;
        coverage = $3;
        sum_coverage += coverage;
        
        if (coverage > 0) {
            covered_bases++;
        }
        
        # Store coverage values for median calculation
        coverages[total_bases] = coverage;
    }
    END {
        if (total_bases == 0) {
            print "Mean coverage: 0";
            print "Median coverage: 0";
            print "Covered bases: 0";
            print "Coverage breadth: 0";
            exit;
        }
        
        mean_coverage = sum_coverage / total_bases;
        coverage_breadth = (covered_bases / total_bases) * 100;
        
        # Calculate median
        n = asort(coverages);
        if (n % 2 == 1) {
            median_coverage = coverages[(n + 1) / 2];
        } else {
            median_coverage = (coverages[n / 2] + coverages[n / 2 + 1]) / 2;
        }
        
        printf "Mean coverage: %.2f\n", mean_coverage;
        printf "Median coverage: %.2f\n", median_coverage;
        printf "Covered bases: %d\n", covered_bases;
        printf "Coverage breadth: %.2f\n", coverage_breadth;
    }' "$coverage_file"
    
    # Clean up coverage file to save space
    rm -f "$coverage_file"
}

# Function to process individual assemblies
process_individual_assemblies() {
    echo "Processing individual assemblies..."
    
    local individual_dir="$ASSEMBLIES_DIR/individual"
    local output_base="$OUTPUT_DIR/individual"
    
    if [[ ! -d "$individual_dir" ]]; then
        echo "Individual assemblies directory not found: $individual_dir"
        return 1
    fi
    
    # Process each assembler
    for assembler in megahit metaspades; do
        local assembler_dir="$individual_dir/$assembler"
        
        if [[ ! -d "$assembler_dir" ]]; then
            echo "  Assembler directory not found: $assembler_dir"
            continue
        fi
        
        echo "  Processing $assembler individual assemblies..."
        
        # Combine all individual assemblies for this assembler
        local combined_file="$output_base/${assembler}_all_contigs.fasta"
        > "$combined_file"
        
        find "$assembler_dir" -name "*_contigs.fasta" | while read -r contig_file; do
            if [[ -f "$contig_file" ]]; then
                cat "$contig_file" >> "$combined_file"
            fi
        done
        
        # Map reads to combined assembly
        local mapping_output="$output_base/$assembler"
        mkdir -p "$mapping_output"
        
        map_reads_to_assembly "$combined_file" "individual_$assembler" "$mapping_output"
    done
}

# Function to process strategic co-assemblies
process_strategic_coassemblies() {
    echo "Processing strategic co-assemblies..."
    
    local strategic_dir="$ASSEMBLIES_DIR/strategic_coassembly/assemblies"
    local output_base="$OUTPUT_DIR/strategic_coassembly"
    
    if [[ ! -d "$strategic_dir" ]]; then
        echo "Strategic co-assemblies directory not found: $strategic_dir"
        return 1
    fi
    
    # Combine all strategic co-assembly contigs
    local combined_file="$output_base/all_strategic_contigs.fasta"
    > "$combined_file"
    
    find "$strategic_dir" -name "*_contigs.fasta" | while read -r contig_file; do
        if [[ -f "$contig_file" ]]; then
            cat "$contig_file" >> "$combined_file"
        fi
    done
    
    # Map reads to combined strategic co-assembly
    mkdir -p "$output_base"
    map_reads_to_assembly "$combined_file" "strategic_coassembly" "$output_base"
}

# Function to process global co-assembly
process_global_coassembly() {
    echo "Processing global co-assembly..."
    
    local global_file="$ASSEMBLIES_DIR/global_coassembly/global_coassembly_contigs.fasta"
    local output_base="$OUTPUT_DIR/global_coassembly"
    
    if [[ ! -f "$global_file" ]]; then
        echo "Global co-assembly file not found: $global_file"
        return 1
    fi
    
    # Map reads to global co-assembly
    mkdir -p "$output_base"
    map_reads_to_assembly "$global_file" "global_coassembly" "$output_base"
}

# Function to process meta-assembly
process_meta_assembly() {
    echo "Processing meta-assembly..."
    
    local meta_file="$ASSEMBLIES_DIR/meta_assembly/meta_assembly_contigs.fasta"
    local output_base="$OUTPUT_DIR/meta_assembly"
    
    if [[ ! -f "$meta_file" ]]; then
        echo "Meta-assembly file not found: $meta_file"
        return 1
    fi
    
    # Map reads to meta-assembly
    mkdir -p "$output_base"
    map_reads_to_assembly "$meta_file" "meta_assembly" "$output_base"
}

# Function to create combined summary report
create_combined_summary() {
    echo "Creating combined read mapping summary..."
    
    local combined_summary="$OUTPUT_DIR/combined/mapping_comparison_summary.txt"
    local combined_csv="$OUTPUT_DIR/combined/mapping_comparison.csv"
    
    mkdir -p "$OUTPUT_DIR/combined"
    
    # Create text summary
    cat > "$combined_summary" << EOF
Read Mapping Analysis - Assembly Strategy Comparison
====================================================

Analysis Date: $(date)
Minimum Mapping Quality: $MIN_MAPQ
Minimum Contig Length: $MIN_CONTIG_LEN bp

Summary by Assembly Strategy:
EOF
    
    # Initialize CSV
    echo "Strategy,Assembler,Total_Reads,Mapped_Reads,Overall_Mapping_Rate,Mean_Coverage_Avg,Coverage_Breadth_Avg" > "$combined_csv"
    
    # Combine results from all strategies
    for strategy_dir in "$OUTPUT_DIR"/{individual,strategic_coassembly,global_coassembly,meta_assembly}; do
        if [[ -d "$strategy_dir" ]]; then
            strategy_name=$(basename "$strategy_dir")
            
            echo "" >> "$combined_summary"
            echo "$strategy_name:" >> "$combined_summary"
            echo "$(printf '%.0s-' {1..50})" >> "$combined_summary"
            
            # Process individual assemblers or single result
            if [[ "$strategy_name" == "individual" ]]; then
                for assembler in megahit metaspades; do
                    local assembler_dir="$strategy_dir/$assembler"
                    if [[ -d "$assembler_dir" ]]; then
                        local stats_file="$assembler_dir/mapping_statistics.txt"
                        local csv_file="$assembler_dir/mapping_summary.csv"
                        
                        if [[ -f "$stats_file" ]]; then
                            echo "" >> "$combined_summary"
                            echo "  $assembler:" >> "$combined_summary"
                            grep -A 20 "Overall Summary:" "$stats_file" >> "$combined_summary" 2>/dev/null || true
                        fi
                        
                        # Extract summary for CSV
                        if [[ -f "$csv_file" ]]; then
                            # Calculate averages from CSV
                            local avg_stats=$(awk -F',' '
                            NR == 1 { next }
                            {
                                total_reads += $2;
                                mapped_reads += $3;
                                sum_coverage += $6;
                                sum_breadth += $8;
                                count++;
                            }
                            END {
                                if (count > 0) {
                                    mapping_rate = (total_reads > 0) ? (mapped_reads * 100 / total_reads) : 0;
                                    avg_coverage = sum_coverage / count;
                                    avg_breadth = sum_breadth / count;
                                    printf "%d,%d,%.2f,%.2f,%.2f", total_reads, mapped_reads, mapping_rate, avg_coverage, avg_breadth;
                                }
                            }' "$csv_file")
                            
                            if [[ -n "$avg_stats" ]]; then
                                echo "individual_$assembler,$assembler,$avg_stats" >> "$combined_csv"
                            fi
                        fi
                    fi
                done
            else
                local stats_file="$strategy_dir/mapping_statistics.txt"
                local csv_file="$strategy_dir/mapping_summary.csv"
                
                if [[ -f "$stats_file" ]]; then
                    echo "" >> "$combined_summary"
                    grep -A 20 "Overall Summary:" "$stats_file" >> "$combined_summary" 2>/dev/null || true
                fi
                
                # Extract summary for CSV
                if [[ -f "$csv_file" ]]; then
                    local avg_stats=$(awk -F',' '
                    NR == 1 { next }
                    {
                        total_reads += $2;
                        mapped_reads += $3;
                        sum_coverage += $6;
                        sum_breadth += $8;
                        count++;
                    }
                    END {
                        if (count > 0) {
                            mapping_rate = (total_reads > 0) ? (mapped_reads * 100 / total_reads) : 0;
                            avg_coverage = sum_coverage / count;
                            avg_breadth = sum_breadth / count;
                            printf "%d,%d,%.2f,%.2f,%.2f", total_reads, mapped_reads, mapping_rate, avg_coverage, avg_breadth;
                        }
                    }' "$csv_file")
                    
                    if [[ -n "$avg_stats" ]]; then
                        echo "$strategy_name,mixed,$avg_stats" >> "$combined_csv"
                    fi
                fi
            fi
        fi
    done
    
    echo "" >> "$combined_summary"
    echo "Files Generated:" >> "$combined_summary"
    echo "  Combined summary: $combined_summary" >> "$combined_summary"
    echo "  CSV comparison: $combined_csv" >> "$combined_summary"
    echo "  Individual results: $OUTPUT_DIR/{individual,strategic_coassembly,global_coassembly,meta_assembly}/" >> "$combined_summary"
    
    echo "Combined summary saved to: $combined_summary"
    echo "CSV comparison saved to: $combined_csv"
}

# Main execution
echo "Starting read mapping analysis pipeline..."
echo ""

# Check dependencies
check_dependencies

# Validate input directories
if [[ ! -d "$ASSEMBLIES_DIR" ]]; then
    echo "Error: Assemblies directory not found: $ASSEMBLIES_DIR"
    exit 1
fi

if [[ ! -d "$READS_DIR" ]]; then
    echo "Error: Reads directory not found: $READS_DIR"
    exit 1
fi

# Find read files
find_read_files
echo ""

# Process each assembly strategy
process_individual_assemblies
echo ""

process_strategic_coassemblies
echo ""

process_global_coassembly
echo ""

process_meta_assembly
echo ""

# Create combined summary
create_combined_summary

# Clean up temporary files
echo ""
echo "Cleaning up temporary files..."
rm -rf "$OUTPUT_DIR/temp"

echo ""
echo "Read mapping analysis completed!"
echo "================================"
echo ""
echo "Results saved to: $OUTPUT_DIR"
echo ""
echo "Key output files:"
echo "  Combined summary: $OUTPUT_DIR/combined/mapping_comparison_summary.txt"
echo "  CSV comparison: $OUTPUT_DIR/combined/mapping_comparison.csv"
echo "  Individual results: $OUTPUT_DIR/{individual,strategic_coassembly,global_coassembly,meta_assembly}/"
echo ""
echo "Next steps:"
echo "  1. Review mapping statistics and coverage patterns"
echo "  2. Run contig statistics analysis: python scripts/quality_assessment/04_contig_stats.py"
echo "  3. Compare all assembly strategies: python scripts/comparison/01_compare_strategies.py"