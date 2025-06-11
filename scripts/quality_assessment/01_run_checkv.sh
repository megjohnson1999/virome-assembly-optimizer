#!/bin/bash

# CheckV quality assessment for viral contigs across all assembly strategies
# Evaluates completeness, contamination, and overall quality of viral contigs

set -euo pipefail

# Configuration
ASSEMBLIES_DIR=${ASSEMBLIES_DIR:-"results/assemblies"}
OUTPUT_DIR=${OUTPUT_DIR:-"results/quality_assessment/checkv"}
CHECKV_DB=${CHECKV_DB:-"/opt/checkv-db-v1.5"}  # Update path as needed
THREADS=${THREADS:-8}
MIN_CONTIG_LEN=${MIN_CONTIG_LEN:-1000}  # CheckV works better with longer contigs

echo "Starting CheckV quality assessment..."
echo "Assemblies directory: $ASSEMBLIES_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "CheckV database: $CHECKV_DB"
echo "Threads: $THREADS"
echo "Minimum contig length: $MIN_CONTIG_LEN"

# Create output directories
mkdir -p "$OUTPUT_DIR"/{individual,strategic_coassembly,global_coassembly,meta_assembly,combined,logs}

# Function to check dependencies
check_dependencies() {
    echo "Checking dependencies..."
    
    if ! command -v checkv &> /dev/null; then
        echo "Error: CheckV not found. Please install CheckV first."
        echo "Installation: conda install -c bioconda checkv"
        exit 1
    fi
    
    if [[ ! -d "$CHECKV_DB" ]]; then
        echo "Error: CheckV database not found at: $CHECKV_DB"
        echo "Please download CheckV database:"
        echo "  checkv download_database /path/to/checkv-db"
        echo "  Then update CHECKV_DB variable in this script"
        exit 1
    fi
    
    echo "Dependencies check passed"
}

# Function to prepare contigs for CheckV
prepare_contigs() {
    local input_file=$1
    local output_file=$2
    local strategy_name=$3
    
    if [[ ! -f "$input_file" ]]; then
        echo "Warning: Input file not found: $input_file"
        return 1
    fi
    
    echo "  Preparing contigs for $strategy_name..."
    
    # Filter by minimum length and clean headers
    seqtk seq -L "$MIN_CONTIG_LEN" "$input_file" | \
    awk -v strategy="$strategy_name" '
    /^>/ {
        # Clean header and add strategy prefix if not present
        header = $0;
        gsub(/[^A-Za-z0-9_.-]/, "_", header);
        if (index(header, strategy) == 0) {
            print ">" strategy "_" substr(header, 2);
        } else {
            print header;
        }
        next;
    }
    {
        # Convert to uppercase and remove invalid characters
        seq = toupper($0);
        gsub(/[^ATCGN]/, "N", seq);
        print seq;
    }' > "$output_file"
    
    local n_contigs=$(grep -c "^>" "$output_file" 2>/dev/null || echo "0")
    echo "    Prepared $n_contigs contigs for CheckV analysis"
    
    return 0
}

# Function to run CheckV on a set of contigs
run_checkv() {
    local input_file=$1
    local output_dir=$2
    local strategy_name=$3
    local log_file="$OUTPUT_DIR/logs/${strategy_name}_checkv.log"
    
    echo "Running CheckV for $strategy_name..."
    echo "  Input: $input_file"
    echo "  Output: $output_dir"
    
    # Remove existing output directory
    if [[ -d "$output_dir" ]]; then
        rm -rf "$output_dir"
    fi
    
    mkdir -p "$output_dir"
    
    # Check if input file has contigs
    local n_contigs=$(grep -c "^>" "$input_file" 2>/dev/null || echo "0")
    
    if [[ $n_contigs -eq 0 ]]; then
        echo "  Warning: No contigs found in $input_file, skipping CheckV"
        return 1
    fi
    
    echo "  Processing $n_contigs contigs..."
    
    # Run CheckV
    local start_time=$(date +%s)
    
    checkv end_to_end \
        "$input_file" \
        "$output_dir" \
        -t "$THREADS" \
        -d "$CHECKV_DB" \
        2>&1 | tee "$log_file"
    
    local exit_code=${PIPESTATUS[0]}
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    if [[ $exit_code -eq 0 ]]; then
        echo "  CheckV completed successfully in ${duration}s"
        
        # Generate summary statistics
        generate_checkv_summary "$output_dir" "$strategy_name"
        
        return 0
    else
        echo "  Error: CheckV failed with exit code $exit_code"
        echo "  Check log file: $log_file"
        return 1
    fi
}

# Function to generate CheckV summary statistics
generate_checkv_summary() {
    local checkv_dir=$1
    local strategy_name=$2
    local summary_file="$checkv_dir/${strategy_name}_summary.txt"
    local quality_summary="$checkv_dir/quality_summary.tsv"
    
    if [[ ! -f "$quality_summary" ]]; then
        echo "Warning: CheckV quality summary not found: $quality_summary"
        return 1
    fi
    
    echo "Generating CheckV summary for $strategy_name..."
    
    # Parse CheckV results using awk
    awk -F'\t' '
    BEGIN {
        print "CheckV Quality Assessment Summary";
        print "=================================";
        print "";
        print "Strategy: " "'$strategy_name'";
        print "Analysis Date: " strftime("%Y-%m-%d %H:%M:%S");
        print "";
        
        # Initialize counters
        total_contigs = 0;
        complete_genomes = 0;
        high_quality = 0;
        medium_quality = 0;
        low_quality = 0;
        not_determined = 0;
        contaminated = 0;
        
        total_length = 0;
        viral_length = 0;
        host_length = 0;
        
        # Quality categories
        split("", quality_counts);
        split("", length_by_quality);
        
        min_completeness = 100;
        max_completeness = 0;
        sum_completeness = 0;
        n_with_completeness = 0;
    }
    
    NR == 1 { 
        # Header line - store column indices
        for (i = 1; i <= NF; i++) {
            if ($i == "contig_id") contig_col = i;
            else if ($i == "contig_length") length_col = i;
            else if ($i == "provirus") provirus_col = i;
            else if ($i == "proviral_length") proviral_length_col = i;
            else if ($i == "gene_count") gene_col = i;
            else if ($i == "viral_genes") viral_genes_col = i;
            else if ($i == "host_genes") host_genes_col = i;
            else if ($i == "checkv_quality") quality_col = i;
            else if ($i == "completeness") completeness_col = i;
            else if ($i == "contamination") contamination_col = i;
        }
        next;
    }
    
    NR > 1 {
        total_contigs++;
        
        # Get values
        contig_length = (length_col ? $length_col : 0);
        quality = (quality_col ? $quality_col : "Not determined");
        completeness = (completeness_col ? $completeness_col : "");
        contamination = (contamination_col ? $contamination_col : "");
        viral_genes = (viral_genes_col ? $viral_genes_col : 0);
        host_genes = (host_genes_col ? $host_genes_col : 0);
        
        total_length += contig_length;
        
        # Count by quality
        quality_counts[quality]++;
        length_by_quality[quality] += contig_length;
        
        # Track completeness
        if (completeness != "" && completeness != "NA") {
            comp_val = parseFloat(completeness);
            if (comp_val >= 0) {
                sum_completeness += comp_val;
                n_with_completeness++;
                if (comp_val < min_completeness) min_completeness = comp_val;
                if (comp_val > max_completeness) max_completeness = comp_val;
            }
        }
        
        # Count contamination
        if (contamination != "" && contamination != "NA") {
            if (parseFloat(contamination) > 0) {
                contaminated++;
            }
        }
        
        # Categorize quality
        if (quality == "Complete") complete_genomes++;
        else if (quality == "High-quality") high_quality++;
        else if (quality == "Medium-quality") medium_quality++;
        else if (quality == "Low-quality") low_quality++;
        else not_determined++;
    }
    
    END {
        if (total_contigs == 0) {
            print "No contigs found in CheckV results";
            exit;
        }
        
        # Overall statistics
        print "Overall Statistics:";
        print "  Total contigs analyzed: " total_contigs;
        print "  Total length: " total_length " bp";
        print "  Mean contig length: " sprintf("%.0f", total_length / total_contigs) " bp";
        print "";
        
        # Quality distribution
        print "Quality Distribution:";
        print "  Complete genomes: " complete_genomes " (" sprintf("%.1f", complete_genomes/total_contigs*100) "%)";
        print "  High-quality: " high_quality " (" sprintf("%.1f", high_quality/total_contigs*100) "%)";
        print "  Medium-quality: " medium_quality " (" sprintf("%.1f", medium_quality/total_contigs*100) "%)";
        print "  Low-quality: " low_quality " (" sprintf("%.1f", low_quality/total_contigs*100) "%)";
        print "  Not determined: " not_determined " (" sprintf("%.1f", not_determined/total_contigs*100) "%)";
        print "";
        
        # Completeness statistics
        if (n_with_completeness > 0) {
            avg_completeness = sum_completeness / n_with_completeness;
            print "Completeness Statistics:";
            print "  Contigs with completeness data: " n_with_completeness;
            print "  Average completeness: " sprintf("%.1f", avg_completeness) "%";
            print "  Min completeness: " sprintf("%.1f", min_completeness) "%";
            print "  Max completeness: " sprintf("%.1f", max_completeness) "%";
            print "";
        }
        
        # Contamination
        print "Contamination:";
        print "  Contaminated contigs: " contaminated " (" sprintf("%.1f", contaminated/total_contigs*100) "%)";
        print "";
        
        # High-quality summary (Complete + High-quality)
        high_qual_total = complete_genomes + high_quality;
        print "High-Quality Summary (Complete + High-quality):";
        print "  Count: " high_qual_total " (" sprintf("%.1f", high_qual_total/total_contigs*100) "%)";
        if (high_qual_total > 0) {
            hq_length = length_by_quality["Complete"] + length_by_quality["High-quality"];
            print "  Total length: " hq_length " bp (" sprintf("%.1f", hq_length/total_length*100) "%)";
            print "  Average length: " sprintf("%.0f", hq_length / high_qual_total) " bp";
        }
    }
    
    function parseFloat(str) {
        # Simple float parsing
        if (str == "" || str == "NA") return -1;
        return str + 0;
    }
    ' "$quality_summary" > "$summary_file"
    
    echo "Summary saved to: $summary_file"
    
    # Create CSV format for easy analysis
    local csv_file="$checkv_dir/${strategy_name}_summary.csv"
    awk -F'\t' '
    BEGIN {
        print "Strategy,Total_Contigs,Complete,High_Quality,Medium_Quality,Low_Quality,Not_Determined,Total_Length,Contaminated";
    }
    NR == 1 { next; }
    {
        strategy = "'$strategy_name'";
        quality[$quality_col]++;
        total++;
        total_len += $length_col;
        if ($contamination_col && $contamination_col != "NA" && $contamination_col > 0) {
            contam++;
        }
    }
    END {
        printf "%s,%d,%d,%d,%d,%d,%d,%d,%d\n", 
               strategy, total, 
               quality["Complete"], quality["High-quality"], 
               quality["Medium-quality"], quality["Low-quality"], 
               quality["Not determined"], total_len, contam;
    }' "$quality_summary" > "$csv_file"
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
        
        # Prepare and run CheckV
        local prepared_file="$output_base/${assembler}_prepared.fasta"
        local checkv_output="$output_base/$assembler"
        
        if prepare_contigs "$combined_file" "$prepared_file" "INDIVIDUAL_${assembler^^}"; then
            run_checkv "$prepared_file" "$checkv_output" "individual_$assembler"
        fi
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
    
    # Prepare and run CheckV
    local prepared_file="$output_base/strategic_prepared.fasta"
    local checkv_output="$output_base/checkv_results"
    
    if prepare_contigs "$combined_file" "$prepared_file" "STRATEGIC"; then
        run_checkv "$prepared_file" "$checkv_output" "strategic_coassembly"
    fi
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
    
    # Prepare and run CheckV
    local prepared_file="$output_base/global_prepared.fasta"
    local checkv_output="$output_base/checkv_results"
    
    if prepare_contigs "$global_file" "$prepared_file" "GLOBAL"; then
        run_checkv "$prepared_file" "$checkv_output" "global_coassembly"
    fi
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
    
    # Prepare and run CheckV
    local prepared_file="$output_base/meta_prepared.fasta"
    local checkv_output="$output_base/checkv_results"
    
    if prepare_contigs "$meta_file" "$prepared_file" "META"; then
        run_checkv "$prepared_file" "$checkv_output" "meta_assembly"
    fi
}

# Function to create combined summary report
create_combined_summary() {
    echo "Creating combined CheckV summary report..."
    
    local combined_summary="$OUTPUT_DIR/combined/checkv_comparison_summary.txt"
    local combined_csv="$OUTPUT_DIR/combined/checkv_comparison.csv"
    
    mkdir -p "$OUTPUT_DIR/combined"
    
    # Create text summary
    cat > "$combined_summary" << EOF
CheckV Quality Assessment - Assembly Strategy Comparison
========================================================

Analysis Date: $(date)
Minimum Contig Length: $MIN_CONTIG_LEN bp
CheckV Database: $CHECKV_DB

Summary by Assembly Strategy:
EOF
    
    # Initialize CSV
    echo "Strategy,Total_Contigs,Complete,High_Quality,Medium_Quality,Low_Quality,Not_Determined,Total_Length,Contaminated,High_Quality_Percentage" > "$combined_csv"
    
    # Combine results from all strategies
    for strategy_dir in "$OUTPUT_DIR"/{individual,strategic_coassembly,global_coassembly,meta_assembly}; do
        if [[ -d "$strategy_dir" ]]; then
            strategy_name=$(basename "$strategy_dir")
            
            echo "" >> "$combined_summary"
            echo "$strategy_name:" >> "$combined_summary"
            echo "$(printf '%.0s-' {1..50})" >> "$combined_summary"
            
            # Find summary files
            find "$strategy_dir" -name "*_summary.txt" | while read -r summary_file; do
                if [[ -f "$summary_file" ]]; then
                    echo "" >> "$combined_summary"
                    cat "$summary_file" >> "$combined_summary"
                fi
            done
            
            # Add to CSV
            find "$strategy_dir" -name "*_summary.csv" | while read -r csv_file; do
                if [[ -f "$csv_file" ]]; then
                    tail -n +2 "$csv_file" | while IFS=, read -r strategy total complete high medium low not_det length contam; do
                        high_quality_total=$((complete + high))
                        high_quality_pct=$(echo "scale=2; $high_quality_total * 100 / $total" | bc -l 2>/dev/null || echo "0")
                        echo "$strategy,$total,$complete,$high,$medium,$low,$not_det,$length,$contam,$high_quality_pct" >> "$combined_csv"
                    done
                fi
            done
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
echo "Starting CheckV quality assessment pipeline..."
echo ""

# Check dependencies
check_dependencies

# Validate assemblies directory
if [[ ! -d "$ASSEMBLIES_DIR" ]]; then
    echo "Error: Assemblies directory not found: $ASSEMBLIES_DIR"
    echo "Please run assembly scripts first"
    exit 1
fi

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

echo ""
echo "CheckV quality assessment completed!"
echo "====================================="
echo ""
echo "Results saved to: $OUTPUT_DIR"
echo ""
echo "Key output files:"
echo "  Combined summary: $OUTPUT_DIR/combined/checkv_comparison_summary.txt"
echo "  CSV comparison: $OUTPUT_DIR/combined/checkv_comparison.csv"
echo "  Individual results: $OUTPUT_DIR/{individual,strategic_coassembly,global_coassembly,meta_assembly}/"
echo ""
echo "Next steps:"
echo "  1. Review CheckV quality summaries"
echo "  2. Run read mapping analysis: bash scripts/quality_assessment/02_read_mapping.sh"
echo "  3. Compare assembly strategies: python scripts/comparison/01_compare_strategies.py"