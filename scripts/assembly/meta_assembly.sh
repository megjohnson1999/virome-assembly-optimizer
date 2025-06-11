#!/bin/bash

# Meta-assembly script - combines contigs from individual assemblies
# Uses overlap-based merging to create a unified assembly

set -euo pipefail

# Configuration
INDIVIDUAL_DIR=${INDIVIDUAL_DIR:-"results/assemblies/individual"}
OUTPUT_DIR=${OUTPUT_DIR:-"results/assemblies/meta_assembly"}
THREADS=${THREADS:-8}
MIN_CONTIG_LEN=${MIN_CONTIG_LEN:-500}
MIN_OVERLAP=${MIN_OVERLAP:-100}
MIN_IDENTITY=${MIN_IDENTITY:-95}
ASSEMBLER_PRIORITY=${ASSEMBLER_PRIORITY:-"metaspades megahit"}  # Priority order for selecting contigs

echo "Starting meta-assembly pipeline..."
echo "Individual assemblies: $INDIVIDUAL_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "Minimum contig length: $MIN_CONTIG_LEN"
echo "Minimum overlap: $MIN_OVERLAP bp"
echo "Minimum identity: $MIN_IDENTITY%"
echo "Assembler priority: $ASSEMBLER_PRIORITY"

# Create output directories
mkdir -p "$OUTPUT_DIR"/{contigs,merged,stats,logs,temp}

# Function to collect individual assembly contigs
collect_individual_contigs() {
    echo "Collecting contigs from individual assemblies..."
    
    local collected_contigs="$OUTPUT_DIR/temp/all_individual_contigs.fasta"
    local contig_info="$OUTPUT_DIR/temp/contig_info.txt"
    
    # Initialize files
    > "$collected_contigs"
    echo "Contig_ID,Sample,Assembler,Length,File_Path" > "$contig_info"
    
    local total_contigs=0
    local total_samples=0
    
    # Process each assembler in priority order
    for assembler in $ASSEMBLER_PRIORITY; do
        local assembler_dir="$INDIVIDUAL_DIR/$assembler"
        
        if [[ ! -d "$assembler_dir" ]]; then
            echo "  Warning: Assembler directory not found: $assembler_dir"
            continue
        fi
        
        echo "  Processing $assembler assemblies..."
        
        # Find all contig files for this assembler
        find "$assembler_dir" -name "*_contigs.fasta" | sort | while read -r contig_file; do
            if [[ ! -f "$contig_file" ]]; then
                continue
            fi
            
            # Extract sample name from filename
            local filename=$(basename "$contig_file")
            local sample_name=${filename%_contigs.fasta}
            
            echo "    Processing: $sample_name"
            
            # Count contigs in this file
            local n_contigs=$(grep -c "^>" "$contig_file" || echo "0")
            
            if [[ $n_contigs -eq 0 ]]; then
                echo "      Warning: No contigs found in $contig_file"
                continue
            fi
            
            # Filter contigs by length and add to collection
            seqtk seq -L "$MIN_CONTIG_LEN" "$contig_file" | \
            awk -v sample="$sample_name" -v assembler="$assembler" -v file="$contig_file" '
            /^>/ {
                # Update contig ID to include meta-assembly prefix
                contig_id = substr($0, 2);  # Remove >
                new_id = "META_" sample "_" assembler "_" contig_id;
                print ">" new_id;
                
                # Store info for tracking
                current_contig = new_id;
                current_length = 0;
            }
            !/^>/ {
                print $0;
                current_length += length($0);
            }
            END {
                if (current_contig != "" && current_length >= '$MIN_CONTIG_LEN') {
                    print current_contig "," sample "," assembler "," current_length "," file >> "'$contig_info'";
                }
            }' >> "$collected_contigs"
            
            ((total_contigs += n_contigs)) || true
            ((total_samples++)) || true
        done
    done
    
    # Count final collected contigs
    local final_contigs=$(grep -c "^>" "$collected_contigs" || echo "0")
    
    echo "  Collection summary:"
    echo "    Total input contigs: $total_contigs"
    echo "    Contigs after length filtering: $final_contigs"
    echo "    Samples processed: $total_samples"
    
    if [[ $final_contigs -eq 0 ]]; then
        echo "Error: No contigs collected for meta-assembly"
        exit 1
    fi
    
    echo "Collected contigs saved to: $collected_contigs"
    return 0
}

# Function to deduplicate contigs using minimap2 and seqkit
deduplicate_contigs() {
    echo "Deduplicating contigs..."
    
    local input_contigs="$OUTPUT_DIR/temp/all_individual_contigs.fasta"
    local deduplicated_contigs="$OUTPUT_DIR/temp/deduplicated_contigs.fasta"
    local overlap_file="$OUTPUT_DIR/temp/overlaps.paf"
    local log_file="$OUTPUT_DIR/logs/deduplication.log"
    
    echo "  Input: $input_contigs"
    echo "  Output: $deduplicated_contigs"
    
    # Step 1: Self-alignment to find overlaps
    echo "  Running self-alignment with minimap2..."
    
    minimap2 \
        -x asm20 \
        -t "$THREADS" \
        --secondary=no \
        -N 50 \
        "$input_contigs" "$input_contigs" > "$overlap_file" 2>"$log_file"
    
    # Step 2: Identify redundant contigs
    echo "  Identifying redundant contigs..."
    
    python3 << 'EOF' > "$OUTPUT_DIR/temp/redundant_contigs.txt"
import sys

# Read overlap results
overlaps = []
with open("results/assemblies/meta_assembly/temp/overlaps.paf", 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 12:
            continue
        
        query = fields[0]
        target = fields[5]
        
        # Skip self-alignments
        if query == target:
            continue
        
        query_len = int(fields[1])
        target_len = int(fields[6])
        match_len = int(fields[10])
        
        # Calculate coverage and identity
        query_cov = match_len / query_len * 100
        target_cov = match_len / target_len * 100
        
        # Extract identity from CIGAR or use match length approximation
        identity = match_len / max(query_len, target_len) * 100
        
        overlaps.append({
            'query': query,
            'target': target,
            'query_len': query_len,
            'target_len': target_len,
            'match_len': match_len,
            'query_cov': query_cov,
            'target_cov': target_cov,
            'identity': identity
        })

# Find redundant contigs
# A contig is redundant if it's entirely contained in another contig
# with high identity and coverage
redundant = set()
min_identity = 95.0
min_coverage = 90.0

for overlap in overlaps:
    # Check if query is contained in target
    if (overlap['query_cov'] >= min_coverage and 
        overlap['identity'] >= min_identity):
        
        # Prefer longer contigs, or metaspades over megahit
        query_sample = overlap['query'].split('_')[1]
        target_sample = overlap['target'].split('_')[1]
        query_assembler = overlap['query'].split('_')[2]
        target_assembler = overlap['target'].split('_')[2]
        
        # If same sample, prefer longer contig
        if query_sample == target_sample:
            if overlap['query_len'] < overlap['target_len']:
                redundant.add(overlap['query'])
            elif overlap['query_len'] > overlap['target_len']:
                redundant.add(overlap['target'])
            else:
                # Same length, prefer metaspades over megahit
                if query_assembler == 'megahit' and target_assembler == 'metaspades':
                    redundant.add(overlap['query'])
                elif query_assembler == 'metaspades' and target_assembler == 'megahit':
                    redundant.add(overlap['target'])
        else:
            # Different samples, prefer longer contig
            if overlap['query_len'] < overlap['target_len']:
                redundant.add(overlap['query'])

# Write redundant contigs
with open("results/assemblies/meta_assembly/temp/redundant_contigs.txt", 'w') as f:
    for contig in redundant:
        f.write(contig + '\n')

print(f"Identified {len(redundant)} redundant contigs")
EOF
    
    local n_redundant=$(wc -l < "$OUTPUT_DIR/temp/redundant_contigs.txt" 2>/dev/null || echo "0")
    echo "  Found $n_redundant redundant contigs"
    
    # Step 3: Remove redundant contigs
    echo "  Removing redundant contigs..."
    
    if [[ $n_redundant -gt 0 ]]; then
        seqtk subseq -v "$input_contigs" "$OUTPUT_DIR/temp/redundant_contigs.txt" > "$deduplicated_contigs"
    else
        cp "$input_contigs" "$deduplicated_contigs"
    fi
    
    local final_count=$(grep -c "^>" "$deduplicated_contigs" || echo "0")
    local original_count=$(grep -c "^>" "$input_contigs" || echo "0")
    
    echo "  Deduplication complete:"
    echo "    Original contigs: $original_count"
    echo "    Redundant contigs: $n_redundant"
    echo "    Final contigs: $final_count"
    echo "    Reduction: $(echo "scale=1; $n_redundant * 100 / $original_count" | bc -l)%"
}

# Function to perform gap-filling and scaffolding (optional)
improve_assembly() {
    echo "Improving meta-assembly..."
    
    local input_contigs="$OUTPUT_DIR/temp/deduplicated_contigs.fasta"
    local improved_contigs="$OUTPUT_DIR/meta_assembly_contigs.fasta"
    
    # For now, just copy the deduplicated contigs
    # In a more advanced version, you could:
    # 1. Run gap filling tools like GapCloser
    # 2. Perform scaffolding
    # 3. Polish with reads
    
    cp "$input_contigs" "$improved_contigs"
    
    echo "Meta-assembly contigs saved to: $improved_contigs"
}

# Function to calculate comprehensive statistics
calculate_meta_assembly_stats() {
    echo "Calculating meta-assembly statistics..."
    
    local contigs_file="$OUTPUT_DIR/meta_assembly_contigs.fasta"
    local stats_file="$OUTPUT_DIR/stats/meta_assembly_stats.txt"
    local comparison_file="$OUTPUT_DIR/stats/comparison_with_individual.txt"
    
    if [[ ! -f "$contigs_file" ]]; then
        echo "Error: Meta-assembly contigs file not found"
        return 1
    fi
    
    # Calculate basic statistics
    awk '
    BEGIN { 
        total_length = 0; 
        contig_count = 0; 
        current_length = 0;
        
        # Track by sample and assembler
        split("", sample_stats);
        split("", assembler_stats);
    }
    /^>/ { 
        if (current_length > 0) {
            lengths[++contig_count] = current_length;
            total_length += current_length;
            
            # Parse contig ID: META_sample_assembler_originalID
            split(prev_header, parts, "_");
            if (length(parts) >= 3) {
                sample = parts[2];
                assembler = parts[3];
                
                sample_stats[sample "_count"]++;
                sample_stats[sample "_length"] += current_length;
                
                assembler_stats[assembler "_count"]++;
                assembler_stats[assembler "_length"] += current_length;
            }
        }
        current_length = 0;
        prev_header = substr($0, 2);  # Remove >
        next; 
    }
    { 
        current_length += length($0); 
    }
    END {
        if (current_length > 0) {
            lengths[++contig_count] = current_length;
            total_length += current_length;
            
            split(prev_header, parts, "_");
            if (length(parts) >= 3) {
                sample = parts[2];
                assembler = parts[3];
                
                sample_stats[sample "_count"]++;
                sample_stats[sample "_length"] += current_length;
                
                assembler_stats[assembler "_count"]++;
                assembler_stats[assembler "_length"] += current_length;
            }
        }
        
        if (contig_count == 0) {
            print "ERROR: No contigs found";
            exit 1;
        }
        
        # Sort lengths
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
        
        # Print main statistics
        print "META-ASSEMBLY STATISTICS";
        print "========================";
        print "";
        print "Overall Assembly:";
        print "  Total Contigs: " contig_count;
        print "  Total Length: " total_length " bp";
        print "  Mean Length: " sprintf("%.2f", total_length / contig_count) " bp";
        print "  Max Length: " sorted_lengths[1] " bp";
        print "  Min Length: " sorted_lengths[n] " bp";
        print "  N50: " n50 " bp";
        print "  N90: " n90 " bp";
        print "";
        
        # Print assembler contribution
        print "Assembler Contributions:";
        for (key in assembler_stats) {
            if (index(key, "_count") > 0) {
                assembler = substr(key, 1, index(key, "_count") - 1);
                count = assembler_stats[key];
                length_key = assembler "_length";
                total_len = assembler_stats[length_key];
                
                print "  " assembler ":";
                print "    Contigs: " count " (" sprintf("%.1f", count/contig_count*100) "%)";
                print "    Length: " total_len " bp (" sprintf("%.1f", total_len/total_length*100) "%)";
                print "    Mean: " sprintf("%.2f", total_len/count) " bp";
            }
        }
        
        # Length distribution
        print "";
        print "Length Distribution:";
        very_long = 0; long_ct = 0; medium = 0; short_ct = 0;
        for (i = 1; i <= n; i++) {
            if (sorted_lengths[i] >= 100000) very_long++;
            else if (sorted_lengths[i] >= 10000) long_ct++;
            else if (sorted_lengths[i] >= 1000) medium++;
            else short_ct++;
        }
        
        print "  Very Long (≥100kb): " very_long " (" sprintf("%.1f", very_long/contig_count*100) "%)";
        print "  Long (10-100kb): " long_ct " (" sprintf("%.1f", long_ct/contig_count*100) "%)";
        print "  Medium (1-10kb): " medium " (" sprintf("%.1f", medium/contig_count*100) "%)";
        print "  Short (<1kb): " short_ct " (" sprintf("%.1f", short_ct/contig_count*100) "%)";
    }' "$contigs_file" > "$stats_file"
    
    echo "Statistics saved to: $stats_file"
    
    # Create CSV summary
    local csv_file="$OUTPUT_DIR/stats/meta_assembly_summary.csv"
    echo "Metric,Value" > "$csv_file"
    grep -E "(Total Contigs|Total Length|Mean Length|Max Length|N50|N90)" "$stats_file" | \
    sed 's/^  //; s/: /,/; s/ bp//; s/ contigs//' >> "$csv_file"
    
    echo "CSV summary saved to: $csv_file"
    
    # Compare with individual assemblies if available
    create_assembly_comparison
}

# Function to create comparison with individual assemblies
create_assembly_comparison() {
    echo "Creating comparison with individual assemblies..."
    
    local comparison_file="$OUTPUT_DIR/stats/individual_vs_meta_comparison.txt"
    
    cat > "$comparison_file" << EOF
Individual vs Meta-Assembly Comparison
======================================

Meta-Assembly Results:
EOF
    
    # Add meta-assembly stats
    grep -A 20 "Overall Assembly:" "$OUTPUT_DIR/stats/meta_assembly_stats.txt" >> "$comparison_file"
    
    cat >> "$comparison_file" << EOF

Individual Assembly Summary:
EOF
    
    # Summarize individual assemblies if available
    if [[ -f "$INDIVIDUAL_DIR/stats/all_samples_summary.txt" ]]; then
        awk '
        NR == 1 { print $0; next }  # Header
        NR > 1 {
            samples++;
            contigs += $3;
            total_length += $4;
            if ($6 > max_n50) max_n50 = $6;
            if (min_n50 == 0 || $6 < min_n50) min_n50 = $6;
            sum_n50 += $6;
        }
        END {
            print "  Total Samples: " samples;
            print "  Total Contigs: " contigs;
            print "  Total Length: " total_length " bp";
            print "  Mean N50: " sprintf("%.0f", sum_n50/samples) " bp";
            print "  Max N50: " max_n50 " bp";
            print "  Min N50: " min_n50 " bp";
        }' "$INDIVIDUAL_DIR/stats/all_samples_summary.csv" >> "$comparison_file"
    else
        echo "  Individual assembly statistics not available" >> "$comparison_file"
    fi
    
    echo "Comparison saved to: $comparison_file"
}

# Function to validate dependencies
check_dependencies() {
    echo "Checking dependencies..."
    
    local missing_deps=()
    
    # Check for required tools
    if ! command -v minimap2 &> /dev/null; then
        missing_deps+=("minimap2")
    fi
    
    if ! command -v seqtk &> /dev/null; then
        missing_deps+=("seqtk")
    fi
    
    if ! command -v python3 &> /dev/null; then
        missing_deps+=("python3")
    fi
    
    if ! command -v bc &> /dev/null; then
        missing_deps+=("bc")
    fi
    
    if [[ ${#missing_deps[@]} -gt 0 ]]; then
        echo "Error: Missing required dependencies:"
        printf '  %s\n' "${missing_deps[@]}"
        echo ""
        echo "Please install missing tools and try again"
        exit 1
    fi
    
    echo "All dependencies found"
}

# Main execution
echo "Starting meta-assembly workflow..."
echo ""

# Check dependencies
check_dependencies

# Validate input directory
if [[ ! -d "$INDIVIDUAL_DIR" ]]; then
    echo "Error: Individual assemblies directory not found: $INDIVIDUAL_DIR"
    echo "Please run individual assemblies first"
    exit 1
fi

# Check if we have any individual assemblies
n_assemblies=0
for assembler in $ASSEMBLER_PRIORITY; do
    if [[ -d "$INDIVIDUAL_DIR/$assembler" ]]; then
        n_assemblies=$((n_assemblies + $(find "$INDIVIDUAL_DIR/$assembler" -name "*_contigs.fasta" | wc -l)))
    fi
done

if [[ $n_assemblies -eq 0 ]]; then
    echo "Error: No individual assembly contigs found in $INDIVIDUAL_DIR"
    echo "Expected directory structure: $INDIVIDUAL_DIR/{megahit,metaspades}/*_contigs.fasta"
    exit 1
fi

echo "Found $n_assemblies individual assemblies to process"
echo ""

# Step 1: Collect contigs from individual assemblies
collect_individual_contigs

# Step 2: Deduplicate overlapping/redundant contigs
deduplicate_contigs

# Step 3: Improve assembly (currently just copies deduplicated contigs)
improve_assembly

# Step 4: Calculate comprehensive statistics
calculate_meta_assembly_stats

# Clean up temporary files
echo ""
echo "Cleaning up temporary files..."
rm -rf "$OUTPUT_DIR/temp"

echo ""
echo "Meta-assembly pipeline completed successfully!"
echo "============================================="
echo ""
echo "Output files:"
echo "  Final contigs: $OUTPUT_DIR/meta_assembly_contigs.fasta"
echo "  Statistics: $OUTPUT_DIR/stats/meta_assembly_stats.txt"
echo "  CSV summary: $OUTPUT_DIR/stats/meta_assembly_summary.csv"
echo "  Comparison: $OUTPUT_DIR/stats/individual_vs_meta_comparison.txt"
echo "  Logs: $OUTPUT_DIR/logs/"
echo ""
echo "Key statistics:"
grep -E "(Total Contigs|Total Length|N50)" "$OUTPUT_DIR/stats/meta_assembly_stats.txt" | sed 's/^  /  /'
echo ""
echo "Next steps:"
echo "  1. Run quality assessment: bash scripts/quality_assessment/01_run_checkv.sh"
echo "  2. Compare with other assembly strategies"
echo "  3. Perform functional annotation"