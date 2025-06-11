#!/usr/bin/env python3

"""
Comprehensive comparison of viral metagenomic assembly strategies.
Integrates results from CheckV, read mapping, and contig statistics to provide
overall strategy recommendations.
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_checkv_results(checkv_dir):
    """Load CheckV quality assessment results."""
    
    logger.info(f"Loading CheckV results from {checkv_dir}")
    
    checkv_results = {}
    
    # Look for combined CSV file
    combined_csv = Path(checkv_dir) / "combined" / "checkv_comparison.csv"
    
    if combined_csv.exists():
        df = pd.read_csv(combined_csv)
        
        # Calculate high-quality percentage
        df['High_Quality_Percentage'] = ((df['Complete'] + df['High_Quality']) / df['Total_Contigs'] * 100).fillna(0)
        
        # Calculate contamination rate
        df['Contamination_Rate'] = (df['Contaminated'] / df['Total_Contigs'] * 100).fillna(0)
        
        for _, row in df.iterrows():
            strategy = row['Strategy']
            checkv_results[strategy] = {
                'total_contigs': row['Total_Contigs'],
                'complete_genomes': row['Complete'],
                'high_quality': row['High_Quality'],
                'medium_quality': row['Medium_Quality'],
                'low_quality': row['Low_Quality'],
                'not_determined': row['Not_Determined'],
                'total_length': row['Total_Length'],
                'contaminated': row['Contaminated'],
                'high_quality_percentage': row['High_Quality_Percentage'],
                'contamination_rate': row['Contamination_Rate']
            }
        
        logger.info(f"Loaded CheckV results for {len(checkv_results)} strategies")
    else:
        logger.warning(f"CheckV combined results not found: {combined_csv}")
    
    return checkv_results

def load_mapping_results(mapping_dir):
    """Load read mapping analysis results."""
    
    logger.info(f"Loading read mapping results from {mapping_dir}")
    
    mapping_results = {}
    
    # Look for combined CSV file
    combined_csv = Path(mapping_dir) / "combined" / "mapping_comparison.csv"
    
    if combined_csv.exists():
        df = pd.read_csv(combined_csv)
        
        for _, row in df.iterrows():
            strategy = row['Strategy']
            mapping_results[strategy] = {
                'total_reads': row['Total_Reads'],
                'mapped_reads': row['Mapped_Reads'],
                'mapping_rate': row['Overall_Mapping_Rate'],
                'mean_coverage_avg': row['Mean_Coverage_Avg'],
                'coverage_breadth_avg': row['Coverage_Breadth_Avg']
            }
        
        logger.info(f"Loaded mapping results for {len(mapping_results)} strategies")
    else:
        logger.warning(f"Mapping combined results not found: {combined_csv}")
    
    return mapping_results

def load_contig_stats(contig_stats_dir):
    """Load contig statistics results."""
    
    logger.info(f"Loading contig statistics from {contig_stats_dir}")
    
    contig_results = {}
    
    # Look for assembly statistics CSV
    stats_csv = Path(contig_stats_dir) / "assembly_statistics.csv"
    
    if stats_csv.exists():
        df = pd.read_csv(stats_csv, index_col=0)
        
        for strategy in df.index:
            row = df.loc[strategy]
            contig_results[strategy] = {
                'n_contigs': row['n_contigs'],
                'total_length': row['total_length'],
                'mean_length': row['mean_length'],
                'n50': row['n50'],
                'n90': row['n90'],
                'l50': row['l50'],
                'gc_content': row['gc_content'],
                'contiguity_ratio': row['contiguity_ratio'],
                'completeness_score': row['completeness_score'],
                'very_long_contigs': row['very_long_contigs'],
                'long_contigs': row['long_contigs']
            }
        
        logger.info(f"Loaded contig statistics for {len(contig_results)} strategies")
    else:
        logger.warning(f"Contig statistics not found: {stats_csv}")
    
    return contig_results

def integrate_results(checkv_results, mapping_results, contig_results):
    """Integrate all quality assessment results into a unified comparison."""
    
    logger.info("Integrating quality assessment results")
    
    # Get all strategies that appear in any dataset
    all_strategies = set()
    all_strategies.update(checkv_results.keys())
    all_strategies.update(mapping_results.keys())
    all_strategies.update(contig_results.keys())
    
    integrated_results = {}
    
    for strategy in all_strategies:
        integrated_results[strategy] = {
            'strategy': strategy,
            # CheckV metrics
            'viral_completeness': 0,
            'high_quality_percentage': 0,
            'contamination_rate': 0,
            'viral_contigs': 0,
            # Mapping metrics
            'mapping_rate': 0,
            'mean_coverage': 0,
            'coverage_breadth': 0,
            # Contig metrics
            'n50': 0,
            'total_length': 0,
            'contiguity_ratio': 0,
            'assembly_efficiency': 0,
            'large_contig_fraction': 0,
            # Overall scores
            'viral_quality_score': 0,
            'assembly_quality_score': 0,
            'mapping_quality_score': 0,
            'overall_score': 0
        }
        
        # Integrate CheckV results
        if strategy in checkv_results:
            cv = checkv_results[strategy]
            integrated_results[strategy].update({
                'viral_completeness': cv['complete_genomes'] + cv['high_quality'],
                'high_quality_percentage': cv['high_quality_percentage'],
                'contamination_rate': cv['contamination_rate'],
                'viral_contigs': cv['total_contigs']
            })
        
        # Integrate mapping results
        if strategy in mapping_results:
            mr = mapping_results[strategy]
            integrated_results[strategy].update({
                'mapping_rate': mr['mapping_rate'],
                'mean_coverage': mr['mean_coverage_avg'],
                'coverage_breadth': mr['coverage_breadth_avg']
            })
        
        # Integrate contig results
        if strategy in contig_results:
            cr = contig_results[strategy]
            large_contigs = cr['very_long_contigs'] + cr['long_contigs']
            large_contig_fraction = (large_contigs / cr['n_contigs'] * 100) if cr['n_contigs'] > 0 else 0
            assembly_efficiency = cr['n_contigs'] / (cr['total_length'] / 1e6) if cr['total_length'] > 0 else 0
            
            integrated_results[strategy].update({
                'n50': cr['n50'],
                'total_length': cr['total_length'],
                'contiguity_ratio': cr['contiguity_ratio'],
                'assembly_efficiency': assembly_efficiency,
                'large_contig_fraction': large_contig_fraction
            })
    
    # Calculate quality scores
    calculate_quality_scores(integrated_results)
    
    return integrated_results

def calculate_quality_scores(integrated_results):
    """Calculate normalized quality scores for comparison."""
    
    logger.info("Calculating quality scores")
    
    strategies = list(integrated_results.keys())
    
    if not strategies:
        return
    
    # Extract values for normalization
    values = {metric: [integrated_results[s][metric] for s in strategies] for s in strategies for metric in integrated_results[s] if isinstance(integrated_results[s][metric], (int, float))}
    
    # Calculate scores for each strategy
    for strategy in strategies:
        result = integrated_results[strategy]
        
        # Viral Quality Score (0-100)
        viral_score = 0
        if result['viral_contigs'] > 0:
            completeness_score = min(result['high_quality_percentage'], 100)
            contamination_penalty = min(result['contamination_rate'] * 2, 50)  # Penalty for contamination
            viral_score = max(0, completeness_score - contamination_penalty)
        
        # Assembly Quality Score (0-100)
        n50_score = min((result['n50'] / 5000) * 50, 50)  # Up to 50 points for N50 >= 5kb
        contiguity_score = min(result['contiguity_ratio'] * 25, 25)  # Up to 25 points for good contiguity
        efficiency_penalty = max(0, (result['assembly_efficiency'] - 2000) / 100)  # Penalty for >2000 contigs/Mbp
        large_contig_score = min(result['large_contig_fraction'], 25)  # Up to 25 points for large contigs
        
        assembly_score = max(0, n50_score + contiguity_score + large_contig_score - efficiency_penalty)
        
        # Mapping Quality Score (0-100)
        mapping_score = (
            min(result['mapping_rate'], 100) * 0.4 +  # Up to 40 points for mapping rate
            min(result['coverage_breadth'], 100) * 0.4 +  # Up to 40 points for coverage breadth
            min(result['mean_coverage'] / 10, 2) * 10  # Up to 20 points for adequate coverage
        )
        
        # Overall Score (weighted average)
        overall_score = (
            viral_score * 0.4 +
            assembly_score * 0.35 +
            mapping_score * 0.25
        )
        
        # Update results
        integrated_results[strategy].update({
            'viral_quality_score': viral_score,
            'assembly_quality_score': assembly_score,
            'mapping_quality_score': mapping_score,
            'overall_score': overall_score
        })

def create_comparison_visualizations(integrated_results, output_dir):
    """Create comprehensive comparison visualizations."""
    
    logger.info("Creating comparison visualizations")
    
    # Convert to DataFrame for easier plotting
    df = pd.DataFrame.from_dict(integrated_results, orient='index')
    
    # Create plots directory
    plots_dir = output_dir / "comparison_plots"
    plots_dir.mkdir(exist_ok=True)
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Figure 1: Overall comparison dashboard
    create_dashboard_plot(df, plots_dir)
    
    # Figure 2: Quality scores comparison
    create_quality_scores_plot(df, plots_dir)
    
    # Figure 3: Detailed metrics comparison
    create_detailed_metrics_plot(df, plots_dir)
    
    # Figure 4: Strategy ranking visualization
    create_ranking_visualization(df, plots_dir)
    
    logger.info(f"Comparison plots saved to {plots_dir}")

def create_dashboard_plot(df, plots_dir):
    """Create overall comparison dashboard."""
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Viral Metagenomic Assembly Strategy Comparison Dashboard', fontsize=16, fontweight='bold')
    
    # Plot 1: Overall scores
    ax1 = axes[0, 0]
    df['overall_score'].plot(kind='bar', ax=ax1, color='steelblue')
    ax1.set_title('Overall Quality Score')
    ax1.set_ylabel('Score (0-100)')
    ax1.tick_params(axis='x', rotation=45)
    
    # Plot 2: Viral quality metrics
    ax2 = axes[0, 1]
    viral_metrics = df[['high_quality_percentage', 'contamination_rate']]
    viral_metrics.plot(kind='bar', ax=ax2)
    ax2.set_title('Viral Quality Metrics')
    ax2.set_ylabel('Percentage')
    ax2.legend(['High Quality %', 'Contamination %'])
    ax2.tick_params(axis='x', rotation=45)
    
    # Plot 3: Assembly metrics
    ax3 = axes[0, 2]
    ax3_twin = ax3.twinx()
    
    df['n50'].plot(kind='bar', ax=ax3, color='orange', alpha=0.7)
    ax3.set_ylabel('N50 (bp)', color='orange')
    ax3.tick_params(axis='y', labelcolor='orange')
    
    (df['total_length'] / 1e6).plot(kind='line', ax=ax3_twin, color='green', marker='o')
    ax3_twin.set_ylabel('Total Length (Mbp)', color='green')
    ax3_twin.tick_params(axis='y', labelcolor='green')
    
    ax3.set_title('Assembly Size Metrics')
    ax3.tick_params(axis='x', rotation=45)
    
    # Plot 4: Mapping performance
    ax4 = axes[1, 0]
    mapping_metrics = df[['mapping_rate', 'coverage_breadth']]
    mapping_metrics.plot(kind='bar', ax=ax4)
    ax4.set_title('Mapping Performance')
    ax4.set_ylabel('Percentage')
    ax4.legend(['Mapping Rate %', 'Coverage Breadth %'])
    ax4.tick_params(axis='x', rotation=45)
    
    # Plot 5: Quality scores breakdown
    ax5 = axes[1, 1]
    quality_scores = df[['viral_quality_score', 'assembly_quality_score', 'mapping_quality_score']]
    quality_scores.plot(kind='bar', stacked=False, ax=ax5)
    ax5.set_title('Quality Scores Breakdown')
    ax5.set_ylabel('Score (0-100)')
    ax5.legend(['Viral', 'Assembly', 'Mapping'])
    ax5.tick_params(axis='x', rotation=45)
    
    # Plot 6: Efficiency vs Quality scatter
    ax6 = axes[1, 2]
    scatter = ax6.scatter(df['assembly_efficiency'], df['overall_score'], 
                         s=df['total_length']/1e5, alpha=0.7, c=df['mapping_rate'], cmap='viridis')
    
    for i, strategy in enumerate(df.index):
        ax6.annotate(strategy, (df.loc[strategy, 'assembly_efficiency'], 
                               df.loc[strategy, 'overall_score']),
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    ax6.set_xlabel('Assembly Efficiency (contigs/Mbp)')
    ax6.set_ylabel('Overall Quality Score')
    ax6.set_title('Efficiency vs Quality\n(size = total length, color = mapping rate)')
    
    plt.colorbar(scatter, ax=ax6, label='Mapping Rate %')
    
    plt.tight_layout()
    plt.savefig(plots_dir / "comparison_dashboard.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_quality_scores_plot(df, plots_dir):
    """Create detailed quality scores comparison."""
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Quality scores radar chart
    ax1 = axes[0, 0]
    create_quality_radar_chart(df, ax1)
    
    # Plot 2: Score components
    ax2 = axes[0, 1]
    score_components = df[['viral_quality_score', 'assembly_quality_score', 'mapping_quality_score']]
    score_components.plot(kind='bar', ax=ax2, stacked=True)
    ax2.set_title('Quality Score Components (Stacked)')
    ax2.set_ylabel('Score')
    ax2.legend(['Viral', 'Assembly', 'Mapping'])
    ax2.tick_params(axis='x', rotation=45)
    
    # Plot 3: Performance vs Score correlation
    ax3 = axes[1, 0]
    metrics = ['high_quality_percentage', 'mapping_rate', 'n50', 'overall_score']
    correlation_data = df[metrics].corr()
    sns.heatmap(correlation_data, annot=True, cmap='coolwarm', center=0, ax=ax3)
    ax3.set_title('Metric Correlations')
    
    # Plot 4: Strategy ranking
    ax4 = axes[1, 1]
    df_sorted = df.sort_values('overall_score', ascending=True)
    colors = plt.cm.RdYlGn(df_sorted['overall_score'] / 100)
    
    bars = ax4.barh(range(len(df_sorted)), df_sorted['overall_score'], color=colors)
    ax4.set_yticks(range(len(df_sorted)))
    ax4.set_yticklabels(df_sorted.index)
    ax4.set_xlabel('Overall Quality Score')
    ax4.set_title('Strategy Ranking')
    
    # Add score labels
    for i, (strategy, score) in enumerate(zip(df_sorted.index, df_sorted['overall_score'])):
        ax4.text(score + 1, i, f'{score:.1f}', va='center')
    
    plt.tight_layout()
    plt.savefig(plots_dir / "quality_scores_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_quality_radar_chart(df, ax):
    """Create radar chart for quality scores comparison."""
    
    # Select top 4 strategies by overall score
    top_strategies = df.nlargest(4, 'overall_score')
    
    if len(top_strategies) == 0:
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center', transform=ax.transAxes)
        return
    
    # Metrics for radar chart
    metrics = ['viral_quality_score', 'assembly_quality_score', 'mapping_quality_score']
    metric_labels = ['Viral Quality', 'Assembly Quality', 'Mapping Quality']
    
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
    colors = ['red', 'blue', 'green', 'orange']
    for i, (strategy, row) in enumerate(top_strategies.iterrows()):
        values = [row[metric] for metric in metrics]
        values += values[:1]  # Complete the circle
        
        ax.plot(angles, values, 'o-', linewidth=2, label=strategy, color=colors[i % len(colors)])
        ax.fill(angles, values, alpha=0.1, color=colors[i % len(colors)])
    
    ax.set_ylim(0, 100)
    ax.set_title('Quality Scores Comparison (Top 4 Strategies)')
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))

def create_detailed_metrics_plot(df, plots_dir):
    """Create detailed metrics comparison plots."""
    
    fig, axes = plt.subplots(3, 2, figsize=(15, 18))
    
    # Plot 1: Viral genome quality
    ax1 = axes[0, 0]
    viral_data = df[['high_quality_percentage', 'contamination_rate']].copy()
    viral_data['clean_percentage'] = 100 - viral_data['contamination_rate']
    
    x = range(len(viral_data))
    width = 0.35
    
    ax1.bar([i - width/2 for i in x], viral_data['high_quality_percentage'], 
           width, label='High Quality %', color='green', alpha=0.7)
    ax1.bar([i + width/2 for i in x], viral_data['clean_percentage'], 
           width, label='Clean (non-contaminated) %', color='blue', alpha=0.7)
    
    ax1.set_xlabel('Strategy')
    ax1.set_ylabel('Percentage')
    ax1.set_title('Viral Genome Quality')
    ax1.set_xticks(x)
    ax1.set_xticklabels(viral_data.index, rotation=45)
    ax1.legend()
    
    # Plot 2: Assembly contiguity
    ax2 = axes[0, 1]
    df['n50'].plot(kind='bar', ax=ax2, color='orange')
    ax2.set_title('Assembly Contiguity (N50)')
    ax2.set_ylabel('N50 (bp)')
    ax2.tick_params(axis='x', rotation=45)
    
    # Plot 3: Coverage analysis
    ax3 = axes[1, 0]
    coverage_data = df[['mapping_rate', 'coverage_breadth', 'mean_coverage']].copy()
    coverage_data['mean_coverage_normalized'] = np.minimum(coverage_data['mean_coverage'] * 10, 100)  # Scale to 0-100
    
    coverage_data[['mapping_rate', 'coverage_breadth', 'mean_coverage_normalized']].plot(kind='bar', ax=ax3)
    ax3.set_title('Coverage Analysis')
    ax3.set_ylabel('Percentage / Normalized Coverage')
    ax3.legend(['Mapping Rate %', 'Coverage Breadth %', 'Mean Coverage (scaled)'])
    ax3.tick_params(axis='x', rotation=45)
    
    # Plot 4: Assembly size comparison
    ax4 = axes[1, 1]
    ax4_twin = ax4.twinx()
    
    (df['total_length'] / 1e6).plot(kind='bar', ax=ax4, color='lightblue', alpha=0.7)
    ax4.set_ylabel('Total Length (Mbp)', color='blue')
    
    df['viral_contigs'].plot(kind='line', ax=ax4_twin, color='red', marker='o', linewidth=2)
    ax4_twin.set_ylabel('Number of Viral Contigs', color='red')
    
    ax4.set_title('Assembly Size vs Viral Content')
    ax4.tick_params(axis='x', rotation=45)
    ax4.tick_params(axis='y', labelcolor='blue')
    ax4_twin.tick_params(axis='y', labelcolor='red')
    
    # Plot 5: Efficiency metrics
    ax5 = axes[2, 0]
    efficiency_data = df[['assembly_efficiency', 'large_contig_fraction']].copy()
    
    x = range(len(efficiency_data))
    width = 0.35
    
    ax5.bar([i - width/2 for i in x], efficiency_data['assembly_efficiency'], 
           width, label='Contigs per Mbp', color='purple', alpha=0.7)
    
    ax5_twin = ax5.twinx()
    ax5_twin.bar([i + width/2 for i in x], efficiency_data['large_contig_fraction'], 
                width, label='Large Contigs %', color='gold', alpha=0.7)
    
    ax5.set_xlabel('Strategy')
    ax5.set_ylabel('Contigs per Mbp', color='purple')
    ax5_twin.set_ylabel('Large Contigs %', color='gold')
    ax5.set_title('Assembly Efficiency')
    ax5.set_xticks(x)
    ax5.set_xticklabels(efficiency_data.index, rotation=45)
    
    # Plot 6: Performance summary heatmap
    ax6 = axes[2, 1]
    
    # Normalize metrics to 0-1 scale for heatmap
    heatmap_metrics = ['overall_score', 'viral_quality_score', 'assembly_quality_score', 
                      'mapping_quality_score', 'high_quality_percentage', 'mapping_rate']
    
    heatmap_data = df[heatmap_metrics].copy()
    for col in heatmap_data.columns:
        max_val = heatmap_data[col].max()
        if max_val > 0:
            heatmap_data[col] = heatmap_data[col] / max_val
    
    sns.heatmap(heatmap_data.T, annot=True, fmt='.2f', cmap='RdYlGn', 
               ax=ax6, cbar_kws={'label': 'Normalized Score'})
    ax6.set_title('Performance Summary Heatmap')
    ax6.set_xlabel('Strategy')
    
    plt.tight_layout()
    plt.savefig(plots_dir / "detailed_metrics_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_ranking_visualization(df, plots_dir):
    """Create comprehensive strategy ranking visualization."""
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Sort by overall score
    df_sorted = df.sort_values('overall_score', ascending=False)
    
    # Plot 1: Overall ranking with score breakdown
    ax1 = axes[0, 0]
    
    score_components = df_sorted[['viral_quality_score', 'assembly_quality_score', 'mapping_quality_score']]
    score_components.plot(kind='barh', stacked=True, ax=ax1)
    
    ax1.set_title('Strategy Ranking with Score Breakdown')
    ax1.set_xlabel('Quality Score')
    ax1.legend(['Viral', 'Assembly', 'Mapping'])
    
    # Add total scores as text
    for i, (strategy, score) in enumerate(zip(df_sorted.index, df_sorted['overall_score'])):
        ax1.text(score + 2, i, f'{score:.1f}', va='center', fontweight='bold')
    
    # Plot 2: Strengths and weaknesses
    ax2 = axes[0, 1]
    
    # Calculate relative strengths (z-scores)
    key_metrics = ['high_quality_percentage', 'mapping_rate', 'n50', 'coverage_breadth']
    strengths_data = df[key_metrics].copy()
    
    # Normalize to z-scores
    for col in strengths_data.columns:
        mean_val = strengths_data[col].mean()
        std_val = strengths_data[col].std()
        if std_val > 0:
            strengths_data[col] = (strengths_data[col] - mean_val) / std_val
    
    sns.heatmap(strengths_data, annot=True, fmt='.1f', cmap='RdBu_r', center=0, ax=ax2)
    ax2.set_title('Relative Strengths/Weaknesses\n(Z-scores)')
    ax2.set_xlabel('Metrics')
    
    # Plot 3: Recommendation confidence
    ax3 = axes[1, 0]
    
    # Calculate confidence based on score gaps and consistency
    scores = df_sorted['overall_score'].values
    
    confidence_scores = []
    for i, score in enumerate(scores):
        if i == 0:  # Best strategy
            gap_to_next = scores[0] - scores[1] if len(scores) > 1 else 20
            confidence = min(100, 50 + gap_to_next * 2)
        else:
            gap_to_prev = scores[i-1] - score
            gap_to_next = score - scores[i+1] if i < len(scores)-1 else 10
            confidence = max(0, 50 - gap_to_prev + gap_to_next)
        
        confidence_scores.append(confidence)
    
    colors = plt.cm.RdYlGn([c/100 for c in confidence_scores])
    bars = ax3.barh(range(len(df_sorted)), confidence_scores, color=colors)
    
    ax3.set_yticks(range(len(df_sorted)))
    ax3.set_yticklabels(df_sorted.index)
    ax3.set_xlabel('Recommendation Confidence (%)')
    ax3.set_title('Strategy Recommendation Confidence')
    
    # Add confidence labels
    for i, conf in enumerate(confidence_scores):
        ax3.text(conf + 2, i, f'{conf:.0f}%', va='center')
    
    # Plot 4: Decision matrix
    ax4 = axes[1, 1]
    
    # Create decision matrix based on use cases
    use_cases = {
        'High Throughput': [0.2, 0.5, 0.3],  # [viral, assembly, mapping] weights
        'Quality Focus': [0.6, 0.3, 0.1],
        'Balanced': [0.4, 0.35, 0.25],
        'Discovery': [0.3, 0.4, 0.3]
    }
    
    decision_scores = {}
    for use_case, weights in use_cases.items():
        scores = []
        for strategy in df.index:
            score = (df.loc[strategy, 'viral_quality_score'] * weights[0] +
                    df.loc[strategy, 'assembly_quality_score'] * weights[1] +
                    df.loc[strategy, 'mapping_quality_score'] * weights[2])
            scores.append(score)
        decision_scores[use_case] = scores
    
    decision_df = pd.DataFrame(decision_scores, index=df.index)
    sns.heatmap(decision_df, annot=True, fmt='.1f', cmap='RdYlGn', ax=ax4)
    ax4.set_title('Use Case Specific Recommendations')
    ax4.set_xlabel('Use Case')
    
    plt.tight_layout()
    plt.savefig(plots_dir / "strategy_ranking_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()

def generate_final_recommendations(integrated_results, output_dir):
    """Generate final comprehensive recommendations."""
    
    logger.info("Generating final recommendations")
    
    # Convert to DataFrame and sort by overall score
    df = pd.DataFrame.from_dict(integrated_results, orient='index')
    df_sorted = df.sort_values('overall_score', ascending=False)
    
    # Generate recommendations report
    recommendations_file = output_dir / "final_assembly_recommendations.txt"
    
    with open(recommendations_file, 'w') as f:
        f.write("FINAL VIRAL METAGENOMIC ASSEMBLY STRATEGY RECOMMENDATIONS\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("EXECUTIVE SUMMARY\n")
        f.write("-" * 17 + "\n")
        
        best_strategy = df_sorted.index[0]
        best_score = df_sorted.loc[best_strategy, 'overall_score']
        
        f.write(f"RECOMMENDED STRATEGY: {best_strategy}\n")
        f.write(f"Overall Quality Score: {best_score:.1f}/100\n\n")
        
        f.write("Key Performance Metrics:\n")
        f.write(f"• Viral Quality: {df_sorted.loc[best_strategy, 'viral_quality_score']:.1f}/100\n")
        f.write(f"• Assembly Quality: {df_sorted.loc[best_strategy, 'assembly_quality_score']:.1f}/100\n")
        f.write(f"• Mapping Quality: {df_sorted.loc[best_strategy, 'mapping_quality_score']:.1f}/100\n")
        f.write(f"• High-Quality Viral Contigs: {df_sorted.loc[best_strategy, 'high_quality_percentage']:.1f}%\n")
        f.write(f"• Assembly N50: {df_sorted.loc[best_strategy, 'n50']:,} bp\n")
        f.write(f"• Read Mapping Rate: {df_sorted.loc[best_strategy, 'mapping_rate']:.1f}%\n\n")
        
        # Strategy ranking
        f.write("COMPLETE STRATEGY RANKING\n")
        f.write("-" * 25 + "\n")
        
        for i, (strategy, row) in enumerate(df_sorted.iterrows(), 1):
            f.write(f"\n{i}. {strategy}\n")
            f.write(f"   Overall Score: {row['overall_score']:.1f}/100\n")
            f.write(f"   Strengths: ")
            
            strengths = []
            if row['viral_quality_score'] >= 70:
                strengths.append("High viral quality")
            if row['assembly_quality_score'] >= 70:
                strengths.append("Good assembly contiguity")
            if row['mapping_quality_score'] >= 70:
                strengths.append("Excellent read mapping")
            if row['n50'] >= 5000:
                strengths.append("Large contigs")
            if row['contamination_rate'] <= 5:
                strengths.append("Low contamination")
            
            f.write(", ".join(strengths) if strengths else "None identified")
            f.write("\n")
            
            f.write(f"   Considerations: ")
            considerations = []
            if row['viral_quality_score'] < 50:
                considerations.append("Low viral completeness")
            if row['assembly_quality_score'] < 50:
                considerations.append("Highly fragmented assembly")
            if row['mapping_quality_score'] < 50:
                considerations.append("Poor read mapping")
            if row['contamination_rate'] > 10:
                considerations.append("High contamination")
            
            f.write(", ".join(considerations) if considerations else "None identified")
            f.write("\n")
        
        # Use case specific recommendations
        f.write("\n\nUSE CASE SPECIFIC RECOMMENDATIONS\n")
        f.write("-" * 34 + "\n")
        
        f.write("\n1. HIGH-QUALITY GENOME RECOVERY:\n")
        f.write("   Prioritize: Viral completeness and low contamination\n")
        viral_best = df_sorted.loc[df_sorted['viral_quality_score'].idxmax()]
        f.write(f"   Recommended: {viral_best.name} (Viral Score: {viral_best['viral_quality_score']:.1f})\n")
        
        f.write("\n2. LARGE-SCALE COMPARATIVE STUDIES:\n")
        f.write("   Prioritize: Read mapping efficiency and consistency\n")
        mapping_best = df_sorted.loc[df_sorted['mapping_quality_score'].idxmax()]
        f.write(f"   Recommended: {mapping_best.name} (Mapping Score: {mapping_best['mapping_quality_score']:.1f})\n")
        
        f.write("\n3. NOVEL VIRUS DISCOVERY:\n")
        f.write("   Prioritize: Assembly contiguity and total length\n")
        assembly_best = df_sorted.loc[df_sorted['assembly_quality_score'].idxmax()]
        f.write(f"   Recommended: {assembly_best.name} (Assembly Score: {assembly_best['assembly_quality_score']:.1f})\n")
        
        f.write("\n4. BALANCED ANALYSIS:\n")
        f.write("   Prioritize: Overall performance across all metrics\n")
        f.write(f"   Recommended: {best_strategy} (Overall Score: {best_score:.1f})\n")
        
        # Implementation guidance
        f.write("\n\nIMPLEMENTATION GUIDANCE\n")
        f.write("-" * 23 + "\n")
        
        f.write(f"\nFor the recommended strategy ({best_strategy}):\n")
        
        if 'individual' in best_strategy.lower():
            f.write("\n• Individual Assembly Strategy:\n")
            f.write("  - Process each sample separately\n")
            f.write("  - Good for heterogeneous datasets\n")
            f.write("  - Lower computational requirements\n")
            f.write("  - May miss low-abundance variants\n")
        
        elif 'strategic' in best_strategy.lower():
            f.write("\n• Strategic Co-assembly Strategy:\n")
            f.write("  - Group similar samples based on k-mer analysis\n")
            f.write("  - Balance between individual and global approaches\n")
            f.write("  - Requires careful sample grouping\n")
            f.write("  - Good for mixed datasets\n")
        
        elif 'global' in best_strategy.lower():
            f.write("\n• Global Co-assembly Strategy:\n")
            f.write("  - Combine all samples together\n")
            f.write("  - Maximum sensitivity for rare variants\n")
            f.write("  - High computational requirements\n")
            f.write("  - Risk of strain mixing\n")
        
        elif 'meta' in best_strategy.lower():
            f.write("\n• Meta-assembly Strategy:\n")
            f.write("  - Combine contigs from individual assemblies\n")
            f.write("  - Reduces redundancy\n")
            f.write("  - Requires overlap detection\n")
            f.write("  - Good compromise approach\n")
        
        # Quality control recommendations
        f.write("\n\nQUALITY CONTROL RECOMMENDATIONS\n")
        f.write("-" * 32 + "\n")
        
        f.write("\n1. Pre-assembly QC:\n")
        f.write("   ✓ Remove adapter sequences and low-quality bases\n")
        f.write("   ✓ Filter host contamination\n")
        f.write("   ✓ Remove PCR duplicates\n")
        f.write("   ✓ Assess k-mer similarity between samples\n")
        
        f.write("\n2. Post-assembly QC:\n")
        f.write("   ✓ Run CheckV for viral genome quality assessment\n")
        f.write("   ✓ Perform read mapping to assess coverage\n")
        f.write("   ✓ Calculate assembly statistics (N50, contiguity)\n")
        f.write("   ✓ Check for contamination and host sequences\n")
        
        f.write("\n3. Validation Steps:\n")
        f.write("   ✓ Compare results across multiple assembly strategies\n")
        f.write("   ✓ Validate high-quality contigs with independent methods\n")
        f.write("   ✓ Perform functional annotation to assess completeness\n")
        f.write("   ✓ Use phylogenetic analysis for taxonomy validation\n")
        
        # Final notes
        f.write("\n\nFINAL NOTES\n")
        f.write("-" * 11 + "\n")
        
        score_gap = df_sorted.iloc[0]['overall_score'] - df_sorted.iloc[1]['overall_score'] if len(df_sorted) > 1 else 0
        
        if score_gap > 20:
            f.write("\n• The recommended strategy shows a clear performance advantage.\n")
        elif score_gap > 10:
            f.write("\n• The recommended strategy shows moderate advantage over alternatives.\n")
        else:
            f.write("\n• Multiple strategies show similar performance. Consider:\n")
            f.write("  - Dataset-specific factors\n")
            f.write("  - Computational resource constraints\n")
            f.write("  - Downstream analysis requirements\n")
        
        f.write(f"\n• This analysis considered {len(df)} assembly strategies.\n")
        f.write("• Recommendations are based on integrated quality metrics.\n")
        f.write("• Consider pilot testing before large-scale implementation.\n")
        f.write("• Regular quality assessment is recommended for ongoing projects.\n")
        
        f.write(f"\n\nGenerated by Viral Metagenomic Assembly Toolkit\n")
        f.write(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # Also save detailed results as CSV
    csv_file = output_dir / "integrated_comparison_results.csv"
    df_sorted.to_csv(csv_file)
    
    logger.info(f"Final recommendations saved to {recommendations_file}")
    logger.info(f"Detailed results saved to {csv_file}")
    
    return df_sorted

def main():
    parser = argparse.ArgumentParser(description="Compare viral metagenomic assembly strategies")
    parser.add_argument("--checkv-dir", default="results/quality_assessment/checkv",
                       help="Directory containing CheckV results")
    parser.add_argument("--mapping-dir", default="results/quality_assessment/read_mapping",
                       help="Directory containing read mapping results")
    parser.add_argument("--contig-stats-dir", default="results/quality_assessment/contig_stats",
                       help="Directory containing contig statistics")
    parser.add_argument("--output-dir", default="results/comparison",
                       help="Output directory for comparison results")
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Load all quality assessment results
        checkv_results = load_checkv_results(args.checkv_dir)
        mapping_results = load_mapping_results(args.mapping_dir)
        contig_results = load_contig_stats(args.contig_stats_dir)
        
        # Integrate results
        integrated_results = integrate_results(checkv_results, mapping_results, contig_results)
        
        if not integrated_results:
            raise ValueError("No integrated results available")
        
        # Create comprehensive visualizations
        create_comparison_visualizations(integrated_results, output_dir)
        
        # Generate final recommendations
        final_ranking = generate_final_recommendations(integrated_results, output_dir)
        
        logger.info("Assembly strategy comparison completed successfully!")
        logger.info(f"Results saved to: {output_dir}")
        
        # Print summary
        print("\nAssembly Strategy Comparison Summary:")
        print("=" * 37)
        print(f"Strategies compared: {len(integrated_results)}")
        
        if len(final_ranking) > 0:
            best_strategy = final_ranking.index[0]
            best_score = final_ranking.loc[best_strategy, 'overall_score']
            print(f"Best strategy: {best_strategy} (score: {best_score:.1f}/100)")
            
            print(f"\nTop 3 strategies:")
            for i, (strategy, row) in enumerate(final_ranking.head(3).iterrows(), 1):
                print(f"{i}. {strategy}: {row['overall_score']:.1f}")
        
        print(f"\nDetailed results: {output_dir}")
        print("See final_assembly_recommendations.txt for complete analysis")
        
    except Exception as e:
        logger.error(f"Assembly strategy comparison failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()