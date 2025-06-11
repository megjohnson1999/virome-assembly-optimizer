#!/usr/bin/env python3

"""
Calculate comprehensive contig statistics for all assembly strategies.
Provides detailed analysis of assembly quality metrics.
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from Bio import SeqIO
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def calculate_assembly_stats(fasta_file):
    """Calculate comprehensive statistics for a FASTA assembly file."""
    
    if not Path(fasta_file).exists():
        logger.warning(f"Assembly file not found: {fasta_file}")
        return None
    
    sequences = []
    total_length = 0
    gc_count = 0
    total_bases = 0
    
    try:
        # Read all sequences
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_len = len(record.seq)
            sequences.append(seq_len)
            total_length += seq_len
            
            # Count GC content
            seq_upper = str(record.seq).upper()
            gc_count += seq_upper.count('G') + seq_upper.count('C')
            total_bases += len(seq_upper.replace('N', ''))
    
    except Exception as e:
        logger.error(f"Error reading {fasta_file}: {e}")
        return None
    
    if not sequences:
        logger.warning(f"No sequences found in {fasta_file}")
        return None
    
    # Sort sequences by length (descending)
    sequences.sort(reverse=True)
    n_contigs = len(sequences)
    
    # Calculate basic statistics
    mean_length = total_length / n_contigs
    median_length = np.median(sequences)
    max_length = sequences[0]
    min_length = sequences[-1]
    
    # Calculate N50, N90, L50, L90
    cumulative = 0
    n50 = n90 = 0
    l50 = l90 = 0
    
    target_50 = total_length * 0.5
    target_90 = total_length * 0.9
    
    for i, length in enumerate(sequences):
        cumulative += length
        
        if cumulative >= target_50 and n50 == 0:
            n50 = length
            l50 = i + 1
        
        if cumulative >= target_90 and n90 == 0:
            n90 = length
            l90 = i + 1
            break
    
    # Calculate GC content
    gc_content = (gc_count / total_bases * 100) if total_bases > 0 else 0
    
    # Length distribution (users should define size categories based on target organisms)
    # Default categories are provided but may need adjustment for viral sequences
    very_short = sum(1 for x in sequences if x < 500)
    short = sum(1 for x in sequences if 500 <= x < 1000)
    medium = sum(1 for x in sequences if 1000 <= x < 5000)
    long_contigs = sum(1 for x in sequences if 5000 <= x < 10000)
    very_long = sum(1 for x in sequences if x >= 10000)
    
    # Quality metrics
    # Assembly contiguity (higher N50 relative to total length is better)
    contiguity_ratio = n50 / (total_length / n_contigs) if n_contigs > 0 else 0
    
    # Assembly completeness proxy (based on contig size distribution)
    completeness_score = (very_long * 5 + long_contigs * 3 + medium * 1) / n_contigs if n_contigs > 0 else 0
    
    return {
        'n_contigs': n_contigs,
        'total_length': total_length,
        'mean_length': mean_length,
        'median_length': median_length,
        'max_length': max_length,
        'min_length': min_length,
        'n50': n50,
        'n90': n90,
        'l50': l50,
        'l90': l90,
        'gc_content': gc_content,
        'very_short_contigs': very_short,
        'short_contigs': short,
        'medium_contigs': medium,
        'long_contigs': long_contigs,
        'very_long_contigs': very_long,
        'contiguity_ratio': contiguity_ratio,
        'completeness_score': completeness_score
    }

def find_assembly_files(assemblies_dir):
    """Find all assembly files across different strategies."""
    
    logger.info(f"Searching for assembly files in {assemblies_dir}")
    
    assembly_files = {}
    
    # Define search patterns for different strategies
    strategies = {
        'individual_megahit': {
            'pattern': 'individual/megahit/*_contigs.fasta',
            'combine': True
        },
        'individual_metaspades': {
            'pattern': 'individual/metaspades/*_contigs.fasta',
            'combine': True
        },
        'strategic_coassembly': {
            'pattern': 'strategic_coassembly/assemblies/*/group_*_contigs.fasta',
            'combine': True
        },
        'global_coassembly': {
            'pattern': 'global_coassembly/global_coassembly_contigs.fasta',
            'combine': False
        },
        'meta_assembly': {
            'pattern': 'meta_assembly/meta_assembly_contigs.fasta',
            'combine': False
        }
    }
    
    base_path = Path(assemblies_dir)
    
    for strategy, config in strategies.items():
        pattern = config['pattern']
        combine = config['combine']
        
        files = list(base_path.glob(pattern))
        
        if files:
            if combine:
                # Combine multiple files into one for analysis
                combined_file = create_combined_assembly(files, strategy, base_path)
                if combined_file:
                    assembly_files[strategy] = combined_file
            else:
                # Single file
                assembly_files[strategy] = files[0]
            
            logger.info(f"Found {strategy}: {len(files)} file(s)")
        else:
            logger.warning(f"No files found for {strategy} (pattern: {pattern})")
    
    return assembly_files

def create_combined_assembly(files, strategy, base_dir):
    """Combine multiple assembly files into a single file for analysis."""
    
    temp_dir = base_dir / "temp_combined"
    temp_dir.mkdir(exist_ok=True)
    
    combined_file = temp_dir / f"{strategy}_combined.fasta"
    
    try:
        with open(combined_file, 'w') as outfile:
            for i, file_path in enumerate(files):
                with open(file_path, 'r') as infile:
                    for line in infile:
                        if line.startswith('>'):
                            # Add file identifier to header
                            sample_id = Path(file_path).stem.replace('_contigs', '')
                            outfile.write(f">{strategy}_{sample_id}_{line[1:]}")
                        else:
                            outfile.write(line)
        
        logger.info(f"Combined {len(files)} files for {strategy}")
        return combined_file
    
    except Exception as e:
        logger.error(f"Failed to combine files for {strategy}: {e}")
        return None

def analyze_all_assemblies(assembly_files, output_dir):
    """Analyze all assembly strategies and create comparison."""
    
    logger.info("Analyzing all assembly strategies")
    
    results = {}
    
    for strategy, file_path in assembly_files.items():
        logger.info(f"Analyzing {strategy}...")
        
        stats = calculate_assembly_stats(file_path)
        if stats:
            stats['strategy'] = strategy
            stats['file_path'] = str(file_path)
            results[strategy] = stats
        else:
            logger.warning(f"Failed to analyze {strategy}")
    
    if not results:
        raise ValueError("No assembly statistics calculated")
    
    # Convert to DataFrame
    df = pd.DataFrame.from_dict(results, orient='index')
    
    # Save detailed results
    csv_file = output_dir / "assembly_statistics.csv"
    df.to_csv(csv_file)
    logger.info(f"Assembly statistics saved to {csv_file}")
    
    return df

def create_assembly_visualizations(stats_df, output_dir):
    """Create comprehensive visualizations of assembly statistics."""
    
    logger.info("Creating assembly visualizations")
    
    # Create plots directory
    plots_dir = output_dir / "assembly_plots"
    plots_dir.mkdir(exist_ok=True)
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Figure 1: Overview comparison
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Plot 1: Number of contigs
    ax1 = axes[0, 0]
    stats_df['n_contigs'].plot(kind='bar', ax=ax1, color='skyblue')
    ax1.set_title('Number of Contigs by Strategy')
    ax1.set_ylabel('Number of Contigs')
    ax1.tick_params(axis='x', rotation=45)
    
    # Plot 2: Total assembly length
    ax2 = axes[0, 1]
    (stats_df['total_length'] / 1e6).plot(kind='bar', ax=ax2, color='lightgreen')
    ax2.set_title('Total Assembly Length by Strategy')
    ax2.set_ylabel('Total Length (Mbp)')
    ax2.tick_params(axis='x', rotation=45)
    
    # Plot 3: N50 comparison
    ax3 = axes[0, 2]
    stats_df['n50'].plot(kind='bar', ax=ax3, color='orange')
    ax3.set_title('N50 by Strategy')
    ax3.set_ylabel('N50 (bp)')
    ax3.tick_params(axis='x', rotation=45)
    
    # Plot 4: Mean contig length
    ax4 = axes[1, 0]
    stats_df['mean_length'].plot(kind='bar', ax=ax4, color='pink')
    ax4.set_title('Mean Contig Length by Strategy')
    ax4.set_ylabel('Mean Length (bp)')
    ax4.tick_params(axis='x', rotation=45)
    
    # Plot 5: GC content
    ax5 = axes[1, 1]
    stats_df['gc_content'].plot(kind='bar', ax=ax5, color='lightcoral')
    ax5.set_title('GC Content by Strategy')
    ax5.set_ylabel('GC Content (%)')
    ax5.tick_params(axis='x', rotation=45)
    
    # Plot 6: Contiguity ratio
    ax6 = axes[1, 2]
    stats_df['contiguity_ratio'].plot(kind='bar', ax=ax6, color='lightsalmon')
    ax6.set_title('Assembly Contiguity by Strategy')
    ax6.set_ylabel('Contiguity Ratio')
    ax6.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(plots_dir / "assembly_overview.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Figure 2: Length distribution analysis
    create_length_distribution_plots(stats_df, plots_dir)
    
    # Figure 3: Quality metrics comparison
    create_quality_metrics_plots(stats_df, plots_dir)
    
    logger.info(f"Assembly plots saved to {plots_dir}")

def create_length_distribution_plots(stats_df, plots_dir):
    """Create detailed length distribution plots."""
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Length distribution categories (stacked bar)
    length_cols = ['very_short_contigs', 'short_contigs', 'medium_contigs', 'long_contigs', 'very_long_contigs']
    length_labels = ['<500bp', '500-1kb', '1-5kb', '5-10kb', '>10kb']
    
    ax1 = axes[0, 0]
    bottom = np.zeros(len(stats_df))
    colors = plt.cm.Set3(np.linspace(0, 1, len(length_cols)))
    
    for i, col in enumerate(length_cols):
        ax1.bar(stats_df.index, stats_df[col], bottom=bottom, 
               label=length_labels[i], color=colors[i])
        bottom += stats_df[col]
    
    ax1.set_title('Contig Length Distribution by Strategy')
    ax1.set_ylabel('Number of Contigs')
    ax1.legend()
    ax1.tick_params(axis='x', rotation=45)
    
    # Plot 2: Length distribution percentages
    ax2 = axes[0, 1]
    length_percentages = stats_df[length_cols].div(stats_df['n_contigs'], axis=0) * 100
    
    bottom = np.zeros(len(stats_df))
    for i, col in enumerate(length_cols):
        ax2.bar(stats_df.index, length_percentages[col], bottom=bottom,
               label=length_labels[i], color=colors[i])
        bottom += length_percentages[col]
    
    ax2.set_title('Contig Length Distribution (Percentage)')
    ax2.set_ylabel('Percentage of Contigs')
    ax2.legend()
    ax2.tick_params(axis='x', rotation=45)
    
    # Plot 3: N50 vs Total Length scatter
    ax3 = axes[1, 0]
    ax3.scatter(stats_df['total_length'] / 1e6, stats_df['n50'], 
               s=100, c=range(len(stats_df)), cmap='viridis', alpha=0.7)
    
    for i, strategy in enumerate(stats_df.index):
        ax3.annotate(strategy, 
                    (stats_df.loc[strategy, 'total_length'] / 1e6, 
                     stats_df.loc[strategy, 'n50']),
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    ax3.set_xlabel('Total Assembly Length (Mbp)')
    ax3.set_ylabel('N50 (bp)')
    ax3.set_title('N50 vs Total Assembly Length')
    
    # Plot 4: L50 comparison
    ax4 = axes[1, 1]
    stats_df['l50'].plot(kind='bar', ax=ax4, color='mediumpurple')
    ax4.set_title('L50 (Number of Contigs for N50)')
    ax4.set_ylabel('L50')
    ax4.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(plots_dir / "length_distribution_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_quality_metrics_plots(stats_df, plots_dir):
    """Create quality metrics comparison plots."""
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Composite quality score
    ax1 = axes[0, 0]
    
    # Calculate composite quality score
    # Normalize metrics to 0-100 scale
    normalized_n50 = (stats_df['n50'] / stats_df['n50'].max()) * 100
    normalized_length = (stats_df['total_length'] / stats_df['total_length'].max()) * 100
    normalized_contiguity = (stats_df['contiguity_ratio'] / stats_df['contiguity_ratio'].max()) * 100
    normalized_completeness = (stats_df['completeness_score'] / stats_df['completeness_score'].max()) * 100
    
    quality_score = (normalized_n50 * 0.3 + 
                    normalized_length * 0.2 + 
                    normalized_contiguity * 0.3 + 
                    normalized_completeness * 0.2)
    
    quality_score.plot(kind='bar', ax=ax1, color='gold')
    ax1.set_title('Composite Quality Score by Strategy')
    ax1.set_ylabel('Quality Score (0-100)')
    ax1.tick_params(axis='x', rotation=45)
    
    # Plot 2: Efficiency metrics (contigs per Mbp)
    ax2 = axes[0, 1]
    efficiency = stats_df['n_contigs'] / (stats_df['total_length'] / 1e6)
    efficiency.plot(kind='bar', ax=ax2, color='lightsteelblue')
    ax2.set_title('Assembly Efficiency (Contigs per Mbp)')
    ax2.set_ylabel('Contigs per Mbp')
    ax2.tick_params(axis='x', rotation=45)
    
    # Plot 3: Large contig content
    ax3 = axes[1, 0]
    large_contig_fraction = (stats_df['very_long_contigs'] + stats_df['long_contigs']) / stats_df['n_contigs'] * 100
    large_contig_fraction.plot(kind='bar', ax=ax3, color='seagreen')
    ax3.set_title('Large Contigs (>5kb) Percentage')
    ax3.set_ylabel('Percentage of Contigs')
    ax3.tick_params(axis='x', rotation=45)
    
    # Plot 4: Assembly metrics radar chart (for best strategies)
    ax4 = axes[1, 1]
    
    # Select top 3 strategies by quality score
    top_strategies = quality_score.nlargest(3)
    
    if len(top_strategies) > 0:
        create_radar_chart(stats_df.loc[top_strategies.index], ax4)
    
    plt.tight_layout()
    plt.savefig(plots_dir / "quality_metrics_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_radar_chart(df, ax):
    """Create a radar chart for comparing top strategies."""
    
    # Define metrics for radar chart
    metrics = ['n50', 'total_length', 'contiguity_ratio', 'completeness_score', 'gc_content']
    metric_labels = ['N50', 'Total Length', 'Contiguity', 'Completeness', 'GC Content']
    
    # Normalize metrics to 0-1 scale
    normalized_df = df[metrics].copy()
    for metric in metrics:
        max_val = normalized_df[metric].max()
        min_val = normalized_df[metric].min()
        if max_val > min_val:
            normalized_df[metric] = (normalized_df[metric] - min_val) / (max_val - min_val)
        else:
            normalized_df[metric] = 1.0
    
    # Number of variables
    N = len(metrics)
    
    # Compute angle for each axis
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]  # Complete the circle
    
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    
    # Add labels
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(metric_labels)
    
    # Plot each strategy
    colors = ['red', 'blue', 'green']
    for i, (strategy, row) in enumerate(normalized_df.iterrows()):
        values = row.tolist()
        values += values[:1]  # Complete the circle
        
        ax.plot(angles, values, 'o-', linewidth=2, label=strategy, color=colors[i % len(colors)])
        ax.fill(angles, values, alpha=0.25, color=colors[i % len(colors)])
    
    ax.set_ylim(0, 1)
    ax.set_title('Assembly Quality Comparison (Top Strategies)')
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))

def generate_assembly_recommendations(stats_df, output_dir):
    """Generate recommendations based on assembly statistics."""
    
    logger.info("Generating assembly strategy recommendations")
    
    # Calculate composite scores for ranking
    scores = {}
    
    for strategy in stats_df.index:
        stats = stats_df.loc[strategy]
        
        # Quality components (normalized 0-100)
        n50_score = min((stats['n50'] / 10000) * 100, 100)  # Good N50 is >10kb
        length_score = min((stats['total_length'] / 1e7) * 100, 100)  # Good total length is >10Mb
        contiguity_score = min(stats['contiguity_ratio'] * 50, 100)  # Good contiguity ratio is >2
        completeness_score = min(stats['completeness_score'] * 20, 100)  # Scale completeness
        
        # Penalties
        fragmentation_penalty = max(0, (stats['n_contigs'] / (stats['total_length'] / 1e6)) - 1000) / 100  # Penalty for >1000 contigs/Mbp
        small_contig_penalty = (stats['very_short_contigs'] + stats['short_contigs']) / stats['n_contigs'] * 50  # Penalty for small contigs
        
        # Final score
        final_score = (n50_score * 0.3 + 
                      length_score * 0.2 + 
                      contiguity_score * 0.25 + 
                      completeness_score * 0.25 - 
                      fragmentation_penalty - 
                      small_contig_penalty)
        
        scores[strategy] = {
            'final_score': max(0, final_score),
            'n50_score': n50_score,
            'length_score': length_score,
            'contiguity_score': contiguity_score,
            'completeness_score': completeness_score,
            'fragmentation_penalty': fragmentation_penalty,
            'small_contig_penalty': small_contig_penalty
        }
    
    # Rank strategies
    ranked_strategies = sorted(scores.items(), key=lambda x: x[1]['final_score'], reverse=True)
    
    # Generate recommendations report
    recommendations_file = output_dir / "assembly_quality_recommendations.txt"
    
    with open(recommendations_file, 'w') as f:
        f.write("Assembly Quality Recommendations Based on Contig Statistics\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("Strategy Ranking (Best to Worst):\n")
        f.write("-" * 34 + "\n")
        
        for i, (strategy, metrics) in enumerate(ranked_strategies, 1):
            stats = stats_df.loc[strategy]
            
            f.write(f"\n{i}. {strategy}\n")
            f.write(f"   Overall Score: {metrics['final_score']:.1f}/100\n")
            f.write(f"   N50: {stats['n50']:,} bp (score: {metrics['n50_score']:.1f})\n")
            f.write(f"   Total Length: {stats['total_length']/1e6:.1f} Mbp (score: {metrics['length_score']:.1f})\n")
            f.write(f"   Contigs: {stats['n_contigs']:,} (efficiency: {stats['n_contigs']/(stats['total_length']/1e6):.0f} per Mbp)\n")
            f.write(f"   Contiguity: {stats['contiguity_ratio']:.2f} (score: {metrics['contiguity_score']:.1f})\n")
            f.write(f"   Large contigs (>5kb): {(stats['very_long_contigs'] + stats['long_contigs'])/stats['n_contigs']*100:.1f}%\n")
        
        f.write("\nDetailed Analysis:\n")
        f.write("-" * 18 + "\n")
        
        best_strategy = ranked_strategies[0][0]
        best_stats = stats_df.loc[best_strategy]
        
        f.write(f"\nRECOMMENDED STRATEGY: {best_strategy}\n")
        f.write(f"This strategy produced the highest quality assembly with:\n")
        f.write(f"- N50 of {best_stats['n50']:,} bp\n")
        f.write(f"- Total length of {best_stats['total_length']/1e6:.1f} Mbp\n")
        f.write(f"- {best_stats['n_contigs']:,} contigs\n")
        f.write(f"- {(best_stats['very_long_contigs'] + best_stats['long_contigs'])/best_stats['n_contigs']*100:.1f}% large contigs (>5kb)\n")
        
        # Strategy-specific observations
        f.write(f"\nStrategy-Specific Observations:\n")
        f.write("-" * 31 + "\n")
        
        for strategy in stats_df.index:
            stats = stats_df.loc[strategy]
            observations = []
            
            if stats['n50'] < 1000:
                observations.append("Low N50 suggests high fragmentation")
            elif stats['n50'] > 10000:
                observations.append("Excellent N50 indicates good contiguity")
            
            if stats['n_contigs'] / (stats['total_length'] / 1e6) > 5000:
                observations.append("High fragmentation (many small contigs)")
            elif stats['n_contigs'] / (stats['total_length'] / 1e6) < 1000:
                observations.append("Good assembly efficiency")
            
            small_frac = (stats['very_short_contigs'] + stats['short_contigs']) / stats['n_contigs']
            if small_frac > 0.5:
                observations.append(f"High proportion of small contigs ({small_frac*100:.1f}%)")
            
            if stats['gc_content'] < 35 or stats['gc_content'] > 65:
                observations.append(f"Unusual GC content ({stats['gc_content']:.1f}%)")
            
            if observations:
                f.write(f"\n{strategy}:\n")
                for obs in observations:
                    f.write(f"  - {obs}\n")
        
        f.write(f"\nGeneral Recommendations:\n")
        f.write("-" * 24 + "\n")
        f.write(f"1. Use {best_strategy} for downstream analyses\n")
        f.write(f"2. Consider the following quality thresholds:\n")
        f.write(f"   - N50 > 5,000 bp (excellent > 10,000 bp)\n")
        f.write(f"   - <2,000 contigs per Mbp\n")
        f.write(f"   - >20% of contigs should be >5kb\n")
        f.write(f"3. Validate assembly quality with:\n")
        f.write(f"   - CheckV for viral genome completeness\n")
        f.write(f"   - Read mapping for coverage uniformity\n")
        f.write(f"   - Functional annotation completeness\n")
    
    logger.info(f"Recommendations saved to {recommendations_file}")
    
    return ranked_strategies

def cleanup_temp_files(assemblies_dir):
    """Clean up temporary combined files."""
    
    temp_dir = Path(assemblies_dir) / "temp_combined"
    if temp_dir.exists():
        import shutil
        shutil.rmtree(temp_dir)
        logger.info("Cleaned up temporary files")

def main():
    parser = argparse.ArgumentParser(description="Calculate comprehensive contig statistics for assembly strategies")
    parser.add_argument("--assemblies-dir", default="results/assemblies",
                       help="Directory containing assembly results")
    parser.add_argument("--output-dir", default="results/quality_assessment/contig_stats",
                       help="Output directory for contig statistics")
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Find assembly files
        assembly_files = find_assembly_files(args.assemblies_dir)
        
        if not assembly_files:
            raise ValueError("No assembly files found")
        
        # Analyze all assemblies
        stats_df = analyze_all_assemblies(assembly_files, output_dir)
        
        # Create visualizations
        create_assembly_visualizations(stats_df, output_dir)
        
        # Generate recommendations
        recommendations = generate_assembly_recommendations(stats_df, output_dir)
        
        # Clean up temporary files
        cleanup_temp_files(args.assemblies_dir)
        
        logger.info("Contig statistics analysis completed successfully!")
        logger.info(f"Results saved to: {output_dir}")
        
        # Print summary
        print("\nContig Statistics Analysis Summary:")
        print("=" * 35)
        print(f"Strategies analyzed: {len(stats_df)}")
        print(f"Best strategy: {recommendations[0][0]} (score: {recommendations[0][1]['final_score']:.1f})")
        print(f"\nResults saved to: {output_dir}")
        
        # Show top 3 strategies
        print("\nTop 3 Strategies:")
        for i, (strategy, score_info) in enumerate(recommendations[:3], 1):
            stats = stats_df.loc[strategy]
            print(f"{i}. {strategy}")
            print(f"   Score: {score_info['final_score']:.1f}")
            print(f"   N50: {stats['n50']:,} bp")
            print(f"   Contigs: {stats['n_contigs']:,}")
            print(f"   Total: {stats['total_length']/1e6:.1f} Mbp")
        
    except Exception as e:
        logger.error(f"Contig statistics analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()