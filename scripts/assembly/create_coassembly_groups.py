#!/usr/bin/env python3

"""
Create co-assembly groups based on k-mer similarity and important variables
identified through PERMANOVA analysis.
"""

import os
import sys
import json
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from itertools import combinations
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_similarity_data(similarity_dir):
    """Load k-mer similarity analysis results."""
    
    logger.info(f"Loading similarity data from {similarity_dir}")
    
    # Load pairwise similarities
    pairwise_file = Path(similarity_dir) / "pairwise_similarities.csv"
    if not pairwise_file.exists():
        raise FileNotFoundError(f"Pairwise similarities file not found: {pairwise_file}")
    
    pairwise_df = pd.read_csv(pairwise_file)
    logger.info(f"Loaded {len(pairwise_df)} pairwise comparisons")
    
    # Load clustering suggestions if available
    clustering_file = Path(similarity_dir) / "clustering_suggestions.json"
    clustering_suggestions = None
    if clustering_file.exists():
        with open(clustering_file, 'r') as f:
            clustering_suggestions = json.load(f)
        logger.info(f"Loaded {len(clustering_suggestions)} clustering suggestions")
    
    return pairwise_df, clustering_suggestions

def load_variable_analysis(variable_dir):
    """Load variable importance analysis results."""
    
    logger.info(f"Loading variable analysis from {variable_dir}")
    
    # Load PERMANOVA results
    permanova_file = Path(variable_dir) / "permanova_results.csv"
    important_vars = []
    
    if permanova_file.exists():
        permanova_df = pd.read_csv(permanova_file)
        
        # Identify significant variables with meaningful effect sizes
        significant_vars = permanova_df[
            (permanova_df['p_value'] <= 0.05) & 
            (permanova_df['R_squared'] >= 0.1)
        ]
        
        important_vars = significant_vars['Variable'].tolist()
        
        logger.info(f"Found {len(important_vars)} important variables: {important_vars}")
    else:
        logger.warning(f"PERMANOVA results not found: {permanova_file}")
    
    return important_vars

def load_metadata(metadata_file):
    """Load sample metadata."""
    
    logger.info(f"Loading metadata from {metadata_file}")
    
    if not Path(metadata_file).exists():
        logger.warning(f"Metadata file not found: {metadata_file}")
        return None
    
    metadata = pd.read_csv(metadata_file)
    logger.info(f"Loaded metadata for {len(metadata)} samples")
    
    return metadata

def create_kmer_based_groups(pairwise_df, similarity_threshold=0.8, max_group_size=8):
    """Create co-assembly groups based on k-mer similarity."""
    
    logger.info(f"Creating k-mer based groups (threshold: {similarity_threshold})")
    
    # Get high similarity pairs
    high_sim_pairs = pairwise_df[pairwise_df['Similarity'] >= similarity_threshold]
    
    if len(high_sim_pairs) == 0:
        logger.info("No high similarity pairs found, recommending individual assemblies")
        # Get all unique samples
        all_samples = set(pairwise_df['Sample1'].tolist() + pairwise_df['Sample2'].tolist())
        return [{'samples': [sample], 'strategy': 'individual'} for sample in all_samples]
    
    # Build groups using graph-like approach
    groups = []
    used_samples = set()
    
    # Create adjacency list of high similarity connections
    connections = {}
    for _, row in high_sim_pairs.iterrows():
        s1, s2 = row['Sample1'], row['Sample2']
        if s1 not in connections:
            connections[s1] = set()
        if s2 not in connections:
            connections[s2] = set()
        connections[s1].add(s2)
        connections[s2].add(s1)
    
    # Find connected components
    def get_connected_component(start_sample, connections, visited):
        component = {start_sample}
        stack = [start_sample]
        visited.add(start_sample)
        
        while stack:
            current = stack.pop()
            for neighbor in connections.get(current, set()):
                if neighbor not in visited and len(component) < max_group_size:
                    visited.add(neighbor)
                    component.add(neighbor)
                    stack.append(neighbor)
        
        return component
    
    visited = set()
    
    for sample in connections:
        if sample not in visited:
            component = get_connected_component(sample, connections, visited)
            if len(component) > 1:
                groups.append({
                    'samples': list(component),
                    'strategy': 'co-assembly',
                    'basis': 'k-mer similarity'
                })
            else:
                groups.append({
                    'samples': list(component),
                    'strategy': 'individual',
                    'basis': 'no similar samples'
                })
    
    # Add samples that weren't in any high similarity pairs
    all_samples = set(pairwise_df['Sample1'].tolist() + pairwise_df['Sample2'].tolist())
    unconnected_samples = all_samples - visited
    
    for sample in unconnected_samples:
        groups.append({
            'samples': [sample],
            'strategy': 'individual',
            'basis': 'no similar samples'
        })
    
    logger.info(f"Created {len(groups)} k-mer based groups")
    return groups

def refine_groups_with_variables(groups, metadata, important_vars, 
                                similarity_within_var=0.9, max_group_size=8):
    """Refine groups using important variables from PERMANOVA analysis."""
    
    if not important_vars or metadata is None:
        logger.info("No important variables or metadata available, keeping k-mer based groups")
        return groups
    
    logger.info(f"Refining groups using variables: {important_vars}")
    
    refined_groups = []
    
    for group in groups:
        if len(group['samples']) <= 1:
            # Single sample groups don't need refinement
            refined_groups.append(group)
            continue
        
        # Get metadata for samples in this group
        group_metadata = metadata[metadata['Sample'].isin(group['samples'])]
        
        if len(group_metadata) == 0:
            logger.warning(f"No metadata found for samples in group: {group['samples']}")
            refined_groups.append(group)
            continue
        
        # Check if samples in group are similar for important variables
        should_split = False
        split_reason = []
        
        for var in important_vars:
            if var not in group_metadata.columns:
                continue
            
            var_values = group_metadata[var].dropna()
            
            if len(var_values) == 0:
                continue
            
            # For categorical variables, check if all samples have same value
            if group_metadata[var].dtype == 'object' or len(var_values.unique()) <= 10:
                unique_values = var_values.unique()
                if len(unique_values) > 1:
                    should_split = True
                    split_reason.append(f"different {var} values: {unique_values}")
            
            # For continuous variables, check coefficient of variation
            else:
                cv = var_values.std() / var_values.mean() if var_values.mean() != 0 else float('inf')
                if cv > 0.3:  # High variability threshold
                    should_split = True
                    split_reason.append(f"high variability in {var} (CV={cv:.2f})")
        
        if should_split:
            # Split group based on important variables
            logger.info(f"Splitting group due to: {'; '.join(split_reason)}")
            
            # For simplicity, convert multi-sample group to individual assemblies
            # In a more sophisticated approach, you could cluster by variable similarity
            for sample in group['samples']:
                refined_groups.append({
                    'samples': [sample],
                    'strategy': 'individual',
                    'basis': f"variable-based split: {'; '.join(split_reason)}"
                })
        else:
            # Keep the group but update the basis
            group['basis'] = f"k-mer similarity + consistent {', '.join(important_vars)}"
            refined_groups.append(group)
    
    logger.info(f"Refined to {len(refined_groups)} groups")
    return refined_groups

def optimize_group_sizes(groups, min_group_size=2, max_group_size=8, target_group_size=4):
    """Optimize group sizes for computational efficiency."""
    
    logger.info("Optimizing group sizes for computational efficiency")
    
    optimized_groups = []
    
    for group in groups:
        if group['strategy'] == 'individual':
            optimized_groups.append(group)
            continue
        
        samples = group['samples']
        n_samples = len(samples)
        
        if n_samples <= max_group_size:
            # Group is already optimal size
            optimized_groups.append(group)
        else:
            # Split large groups
            logger.info(f"Splitting large group of {n_samples} samples")
            
            n_subgroups = (n_samples + target_group_size - 1) // target_group_size
            subgroup_size = n_samples // n_subgroups
            
            for i in range(n_subgroups):
                start_idx = i * subgroup_size
                if i == n_subgroups - 1:
                    # Last subgroup gets remaining samples
                    subgroup_samples = samples[start_idx:]
                else:
                    subgroup_samples = samples[start_idx:start_idx + subgroup_size]
                
                optimized_groups.append({
                    'samples': subgroup_samples,
                    'strategy': 'co-assembly',
                    'basis': group.get('basis', 'k-mer similarity') + ' (size-optimized)'
                })
    
    logger.info(f"Optimized to {len(optimized_groups)} groups")
    return optimized_groups

def estimate_computational_requirements(groups, avg_reads_per_sample=10000000):
    """Estimate computational requirements for each group."""
    
    logger.info("Estimating computational requirements")
    
    for i, group in enumerate(groups):
        n_samples = len(group['samples'])
        
        # Estimate total reads
        estimated_reads = n_samples * avg_reads_per_sample
        
        # Estimate memory requirements (rough estimates)
        if group['strategy'] == 'individual':
            est_memory_gb = min(8, max(4, estimated_reads / 5000000))
            est_time_hours = min(4, max(1, estimated_reads / 10000000))
        else:  # co-assembly
            est_memory_gb = min(64, max(16, estimated_reads / 2000000))
            est_time_hours = min(12, max(2, estimated_reads / 5000000))
        
        group['estimated_memory_gb'] = int(est_memory_gb)
        group['estimated_time_hours'] = int(est_time_hours)
        group['estimated_reads'] = estimated_reads
        
        groups[i] = group
    
    return groups

def generate_assembly_commands(groups, output_dir, input_dir="data/cleaned_reads"):
    """Generate shell commands for each assembly group."""
    
    logger.info("Generating assembly commands")
    
    commands_dir = Path(output_dir) / "commands"
    commands_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate individual commands for each group
    for i, group in enumerate(groups):
        group_id = f"group_{i+1:03d}"
        samples = group['samples']
        strategy = group['strategy']
        
        if strategy == 'individual':
            for j, sample in enumerate(samples):
                cmd_file = commands_dir / f"{group_id}_{sample}_individual.sh"
                
                with open(cmd_file, 'w') as f:
                    f.write("#!/bin/bash\n\n")
                    f.write(f"# Individual assembly for sample: {sample}\n")
                    f.write(f"# Estimated memory: {group['estimated_memory_gb']}GB\n")
                    f.write(f"# Estimated time: {group['estimated_time_hours']} hours\n\n")
                    
                    f.write(f"SAMPLE_NAME='{sample}'\n")
                    f.write(f"INPUT_DIR='{input_dir}'\n")
                    f.write(f"OUTPUT_DIR='results/assemblies/strategic_coassembly/{group_id}'\n")
                    f.write(f"THREADS=8\n")
                    f.write(f"MEMORY={group['estimated_memory_gb']}\n\n")
                    
                    f.write("# Create output directory\n")
                    f.write("mkdir -p $OUTPUT_DIR\n\n")
                    
                    f.write("# Find input files\n")
                    f.write("R1_FILE=$(find $INPUT_DIR -name \"${SAMPLE_NAME}_*R1*.fastq*\" | head -1)\n")
                    f.write("R2_FILE=$(find $INPUT_DIR -name \"${SAMPLE_NAME}_*R2*.fastq*\" | head -1)\n\n")
                    
                    f.write("# Run MEGAHIT\n")
                    f.write("megahit \\\n")
                    f.write("  -1 $R1_FILE \\\n")
                    f.write("  -2 $R2_FILE \\\n")
                    f.write("  -o $OUTPUT_DIR/megahit_${SAMPLE_NAME} \\\n")
                    f.write("  --num-cpu-threads $THREADS \\\n")
                    f.write("  --memory 0.8 \\\n")
                    f.write("  --min-contig-len 500 \\\n")
                    f.write("  --presets meta-sensitive\n\n")
                    
                    f.write("# Copy and rename final contigs\n")
                    f.write("cp $OUTPUT_DIR/megahit_${SAMPLE_NAME}/final.contigs.fa \\\n")
                    f.write("   $OUTPUT_DIR/${SAMPLE_NAME}_contigs.fasta\n\n")
                    
                    f.write("# Add sample prefix to contig names\n")
                    f.write(f"sed -i 's/^>/>{sample}_/' $OUTPUT_DIR/${{SAMPLE_NAME}}_contigs.fasta\n")
                
                cmd_file.chmod(0o755)
        
        else:  # co-assembly
            cmd_file = commands_dir / f"{group_id}_coassembly.sh"
            
            with open(cmd_file, 'w') as f:
                f.write("#!/bin/bash\n\n")
                f.write(f"# Co-assembly for group: {group_id}\n")
                f.write(f"# Samples: {', '.join(samples)}\n")
                f.write(f"# Basis: {group.get('basis', 'k-mer similarity')}\n")
                f.write(f"# Estimated memory: {group['estimated_memory_gb']}GB\n")
                f.write(f"# Estimated time: {group['estimated_time_hours']} hours\n\n")
                
                f.write("# Configuration\n")
                f.write(f"GROUP_ID='{group_id}'\n")
                f.write(f"INPUT_DIR='{input_dir}'\n")
                f.write(f"OUTPUT_DIR='results/assemblies/strategic_coassembly/{group_id}'\n")
                f.write(f"THREADS=8\n")
                f.write(f"MEMORY={group['estimated_memory_gb']}\n\n")
                
                f.write("# Sample list\n")
                f.write("SAMPLES=(\n")
                for sample in samples:
                    f.write(f"  '{sample}'\n")
                f.write(")\n\n")
                
                f.write("# Create output directory\n")
                f.write("mkdir -p $OUTPUT_DIR\n\n")
                
                f.write("# Find and concatenate input files\n")
                f.write("R1_FILES=()\n")
                f.write("R2_FILES=()\n\n")
                
                f.write("for sample in \"${SAMPLES[@]}\"; do\n")
                f.write("  r1_file=$(find $INPUT_DIR -name \"${sample}_*R1*.fastq*\" | head -1)\n")
                f.write("  r2_file=$(find $INPUT_DIR -name \"${sample}_*R2*.fastq*\" | head -1)\n")
                f.write("  \n")
                f.write("  if [[ -f \"$r1_file\" && -f \"$r2_file\" ]]; then\n")
                f.write("    R1_FILES+=(\"$r1_file\")\n")
                f.write("    R2_FILES+=(\"$r2_file\")\n")
                f.write("    echo \"Found files for $sample: $r1_file, $r2_file\"\n")
                f.write("  else\n")
                f.write("    echo \"Warning: Files not found for $sample\"\n")
                f.write("  fi\n")
                f.write("done\n\n")
                
                f.write("# Check if we have files to process\n")
                f.write("if [[ ${#R1_FILES[@]} -eq 0 ]]; then\n")
                f.write("  echo \"Error: No input files found\"\n")
                f.write("  exit 1\n")
                f.write("fi\n\n")
                
                f.write("# Run MEGAHIT co-assembly\n")
                f.write("megahit \\\n")
                f.write("  -1 $(IFS=','; echo \"${R1_FILES[*]}\") \\\n")
                f.write("  -2 $(IFS=','; echo \"${R2_FILES[*]}\") \\\n")
                f.write("  -o $OUTPUT_DIR/megahit_coassembly \\\n")
                f.write("  --num-cpu-threads $THREADS \\\n")
                f.write("  --memory 0.8 \\\n")
                f.write("  --min-contig-len 500 \\\n")
                f.write("  --presets meta-sensitive\n\n")
                
                f.write("# Copy and rename final contigs\n")
                f.write("cp $OUTPUT_DIR/megahit_coassembly/final.contigs.fa \\\n")
                f.write("   $OUTPUT_DIR/${GROUP_ID}_contigs.fasta\n\n")
                
                f.write("# Add group prefix to contig names\n")
                f.write(f"sed -i 's/^>/>{group_id}_/' $OUTPUT_DIR/${{GROUP_ID}}_contigs.fasta\n")
            
            cmd_file.chmod(0o755)
    
    # Generate master script to run all assemblies
    master_script = commands_dir / "run_all_assemblies.sh"
    with open(master_script, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write("# Master script to run all strategic co-assemblies\n")
        f.write("# Generated automatically by create_coassembly_groups.py\n\n")
        
        f.write("SCRIPT_DIR=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"\n")
        f.write("LOG_DIR=\"results/assemblies/strategic_coassembly/logs\"\n")
        f.write("mkdir -p \"$LOG_DIR\"\n\n")
        
        f.write("echo \"Starting strategic co-assembly pipeline...\"\n")
        f.write("echo \"Total groups: $(ls $SCRIPT_DIR/group_*.sh | wc -l)\"\n\n")
        
        # List all command scripts
        command_scripts = list(commands_dir.glob("group_*.sh"))
        for script in sorted(command_scripts):
            script_name = script.name
            log_name = script_name.replace('.sh', '.log')
            
            f.write(f"echo \"Running {script_name}...\"\n")
            f.write(f"bash \"$SCRIPT_DIR/{script_name}\" 2>&1 | tee \"$LOG_DIR/{log_name}\"\n")
            f.write(f"echo \"Completed {script_name}\"\n\n")
        
        f.write("echo \"All strategic co-assemblies completed!\"\n")
    
    master_script.chmod(0o755)
    
    logger.info(f"Generated {len(list(commands_dir.glob('group_*.sh')))} assembly commands")
    logger.info(f"Master script: {master_script}")

def main():
    parser = argparse.ArgumentParser(description="Create co-assembly groups for viral metagenomics")
    parser.add_argument("--similarity-dir", default="results/similarity_analysis",
                       help="Directory containing k-mer similarity analysis results")
    parser.add_argument("--variable-dir", default="results/variable_analysis",
                       help="Directory containing variable importance analysis results")
    parser.add_argument("--metadata", default="examples/metadata_template.csv",
                       help="Sample metadata file")
    parser.add_argument("--output-dir", default="results/coassembly_groups",
                       help="Output directory for co-assembly group definitions")
    parser.add_argument("--similarity-threshold", type=float, default=0.8,
                       help="K-mer similarity threshold for grouping")
    parser.add_argument("--max-group-size", type=int, default=8,
                       help="Maximum samples per co-assembly group")
    parser.add_argument("--input-dir", default="data/cleaned_reads",
                       help="Directory containing cleaned sequencing reads")
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Load k-mer similarity data
        pairwise_df, clustering_suggestions = load_similarity_data(args.similarity_dir)
        
        # Load variable importance analysis
        important_vars = load_variable_analysis(args.variable_dir)
        
        # Load metadata
        metadata = load_metadata(args.metadata)
        
        # Create initial groups based on k-mer similarity
        groups = create_kmer_based_groups(
            pairwise_df, 
            similarity_threshold=args.similarity_threshold,
            max_group_size=args.max_group_size
        )
        
        # Refine groups using important variables
        groups = refine_groups_with_variables(groups, metadata, important_vars)
        
        # Optimize group sizes
        groups = optimize_group_sizes(groups, max_group_size=args.max_group_size)
        
        # Estimate computational requirements
        groups = estimate_computational_requirements(groups)
        
        # Save group definitions
        groups_file = output_dir / "coassembly_groups.json"
        with open(groups_file, 'w') as f:
            json.dump(groups, f, indent=2)
        
        logger.info(f"Co-assembly groups saved to: {groups_file}")
        
        # Generate summary
        summary = {
            'total_groups': len(groups),
            'individual_assemblies': sum(1 for g in groups if g['strategy'] == 'individual'),
            'co_assemblies': sum(1 for g in groups if g['strategy'] == 'co-assembly'),
            'total_samples': sum(len(g['samples']) for g in groups),
            'largest_group_size': max(len(g['samples']) for g in groups),
            'estimated_total_memory_gb': sum(g['estimated_memory_gb'] for g in groups),
            'estimated_total_time_hours': sum(g['estimated_time_hours'] for g in groups)
        }
        
        summary_file = output_dir / "assembly_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Create human-readable summary
        summary_text_file = output_dir / "assembly_summary.txt"
        with open(summary_text_file, 'w') as f:
            f.write("Strategic Co-assembly Groups Summary\n")
            f.write("====================================\n\n")
            f.write(f"Total groups: {summary['total_groups']}\n")
            f.write(f"Individual assemblies: {summary['individual_assemblies']}\n")
            f.write(f"Co-assemblies: {summary['co_assemblies']}\n")
            f.write(f"Total samples: {summary['total_samples']}\n")
            f.write(f"Largest group size: {summary['largest_group_size']}\n")
            f.write(f"Estimated total memory: {summary['estimated_total_memory_gb']} GB\n")
            f.write(f"Estimated total time: {summary['estimated_total_time_hours']} hours\n\n")
            
            f.write("Group Details:\n")
            f.write("-" * 50 + "\n")
            
            for i, group in enumerate(groups, 1):
                f.write(f"Group {i:03d}: {group['strategy']}\n")
                f.write(f"  Samples ({len(group['samples'])}): {', '.join(group['samples'])}\n")
                f.write(f"  Basis: {group.get('basis', 'not specified')}\n")
                f.write(f"  Memory: {group['estimated_memory_gb']} GB\n")
                f.write(f"  Time: {group['estimated_time_hours']} hours\n\n")
        
        # Generate assembly commands
        generate_assembly_commands(groups, args.output_dir, args.input_dir)
        
        logger.info("Co-assembly group creation completed successfully!")
        logger.info(f"Summary: {summary}")
        
    except Exception as e:
        logger.error(f"Error creating co-assembly groups: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()