#!/usr/bin/env python3

"""
Analyze coverage patterns from read mapping results to assess assembly quality.
Creates detailed coverage statistics and visualizations.
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import subprocess
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_mapping_results(mapping_dir):
    """Load read mapping results from multiple strategies."""
    
    logger.info(f"Loading mapping results from {mapping_dir}")
    
    results = {}
    
    # Define strategy directories
    strategies = {
        'individual_megahit': 'individual/megahit',
        'individual_metaspades': 'individual/metaspades', 
        'strategic_coassembly': 'strategic_coassembly',
        'global_coassembly': 'global_coassembly',
        'meta_assembly': 'meta_assembly'
    }
    
    for strategy, subdir in strategies.items():
        strategy_dir = Path(mapping_dir) / subdir
        csv_file = strategy_dir / "mapping_summary.csv"
        
        if csv_file.exists():
            try:
                df = pd.read_csv(csv_file)
                df['Strategy'] = strategy
                results[strategy] = df
                logger.info(f"Loaded {len(df)} samples for {strategy}")
            except Exception as e:
                logger.warning(f"Failed to load {csv_file}: {e}")
        else:
            logger.warning(f"Mapping results not found: {csv_file}")
    
    if not results:
        raise ValueError("No mapping results found")
    
    # Combine all results
    combined_df = pd.concat(results.values(), ignore_index=True)
    logger.info(f"Combined {len(combined_df)} total mapping results")
    
    return combined_df, results

def calculate_coverage_distributions(mapping_results):
    """Calculate coverage distribution statistics across strategies."""
    
    logger.info("Calculating coverage distribution statistics")
    
    # Group by strategy and calculate statistics
    coverage_stats = mapping_results.groupby('Strategy').agg({
        'Mean_Coverage': ['count', 'mean', 'median', 'std', 'min', 'max'],
        'Median_Coverage': ['mean', 'median', 'std', 'min', 'max'],
        'Coverage_Breadth': ['mean', 'median', 'std', 'min', 'max'],
        'Mapping_Rate': ['mean', 'median', 'std', 'min', 'max']
    }).round(2)
    
    # Flatten column names
    coverage_stats.columns = ['_'.join(col).strip() for col in coverage_stats.columns]
    
    return coverage_stats

def create_coverage_visualizations(mapping_results, output_dir):
    """Create comprehensive coverage visualizations."""
    
    logger.info("Creating coverage visualizations")
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure directory
    plots_dir = Path(output_dir) / "coverage_plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Mapping rate comparison
    plt.figure(figsize=(12, 8))
    
    # Box plot of mapping rates
    plt.subplot(2, 2, 1)
    sns.boxplot(data=mapping_results, x='Strategy', y='Mapping_Rate')
    plt.title('Mapping Rate Distribution by Strategy')
    plt.ylabel('Mapping Rate (%)')
    plt.xticks(rotation=45)
    
    # 2. Mean coverage comparison
    plt.subplot(2, 2, 2)
    sns.boxplot(data=mapping_results, x='Strategy', y='Mean_Coverage')
    plt.title('Mean Coverage Distribution by Strategy')
    plt.ylabel('Mean Coverage (x)')
    plt.xticks(rotation=45)
    plt.yscale('log')
    
    # 3. Coverage breadth comparison
    plt.subplot(2, 2, 3)
    sns.boxplot(data=mapping_results, x='Strategy', y='Coverage_Breadth')
    plt.title('Coverage Breadth Distribution by Strategy')
    plt.ylabel('Coverage Breadth (%)')
    plt.xticks(rotation=45)
    
    # 4. Mapping rate vs coverage scatter
    plt.subplot(2, 2, 4)
    for strategy in mapping_results['Strategy'].unique():
        strategy_data = mapping_results[mapping_results['Strategy'] == strategy]
        plt.scatter(strategy_data['Mapping_Rate'], strategy_data['Mean_Coverage'], 
                   label=strategy, alpha=0.7)
    
    plt.xlabel('Mapping Rate (%)')
    plt.ylabel('Mean Coverage (x)')
    plt.title('Mapping Rate vs Mean Coverage')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.yscale('log')
    
    plt.tight_layout()
    plt.savefig(plots_dir / "coverage_overview.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Individual detailed plots
    create_detailed_coverage_plots(mapping_results, plots_dir)
    
    logger.info(f"Coverage plots saved to {plots_dir}")

def create_detailed_coverage_plots(mapping_results, plots_dir):
    """Create detailed coverage analysis plots."""
    
    # 1. Coverage uniformity analysis
    plt.figure(figsize=(14, 10))
    
    # Calculate coverage uniformity (ratio of median to mean)
    mapping_results['Coverage_Uniformity'] = mapping_results['Median_Coverage'] / mapping_results['Mean_Coverage']
    mapping_results['Coverage_Uniformity'] = mapping_results['Coverage_Uniformity'].fillna(0)
    
    # Plot 1: Coverage uniformity
    plt.subplot(2, 3, 1)
    sns.boxplot(data=mapping_results, x='Strategy', y='Coverage_Uniformity')
    plt.title('Coverage Uniformity (Median/Mean)')
    plt.ylabel('Uniformity Ratio')
    plt.xticks(rotation=45)
    
    # Plot 2: Coverage breadth vs mapping rate
    plt.subplot(2, 3, 2)
    sns.scatterplot(data=mapping_results, x='Mapping_Rate', y='Coverage_Breadth', 
                   hue='Strategy', alpha=0.7)
    plt.title('Coverage Breadth vs Mapping Rate')
    plt.xlabel('Mapping Rate (%)')
    plt.ylabel('Coverage Breadth (%)')
    
    # Plot 3: Mean vs median coverage
    plt.subplot(2, 3, 3)
    sns.scatterplot(data=mapping_results, x='Mean_Coverage', y='Median_Coverage', 
                   hue='Strategy', alpha=0.7)
    plt.title('Mean vs Median Coverage')
    plt.xlabel('Mean Coverage (x)')
    plt.ylabel('Median Coverage (x)')
    plt.xscale('log')
    plt.yscale('log')
    
    # Plot 4: Distribution of mapping rates
    plt.subplot(2, 3, 4)
    for strategy in mapping_results['Strategy'].unique():
        strategy_data = mapping_results[mapping_results['Strategy'] == strategy]
        plt.hist(strategy_data['Mapping_Rate'], alpha=0.7, label=strategy, bins=20)
    plt.xlabel('Mapping Rate (%)')
    plt.ylabel('Frequency')
    plt.title('Distribution of Mapping Rates')
    plt.legend()
    
    # Plot 5: Coverage quality score
    plt.subplot(2, 3, 5)
    # Calculate a composite quality score
    mapping_results['Quality_Score'] = (
        mapping_results['Mapping_Rate'] * 0.4 +
        mapping_results['Coverage_Breadth'] * 0.4 +
        (100 - abs(mapping_results['Coverage_Uniformity'] - 1) * 100) * 0.2
    )
    
    sns.boxplot(data=mapping_results, x='Strategy', y='Quality_Score')
    plt.title('Composite Quality Score')
    plt.ylabel('Quality Score')
    plt.xticks(rotation=45)
    
    # Plot 6: Strategy ranking
    plt.subplot(2, 3, 6)
    strategy_means = mapping_results.groupby('Strategy')['Quality_Score'].mean().sort_values(ascending=False)
    strategy_means.plot(kind='bar')
    plt.title('Average Quality Score by Strategy')
    plt.ylabel('Average Quality Score')
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig(plots_dir / "detailed_coverage_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()

def analyze_coverage_patterns(mapping_results, output_dir):
    """Analyze patterns in coverage data to identify assembly quality issues."""
    
    logger.info("Analyzing coverage patterns")
    
    analysis_results = {}
    
    # 1. Identify samples with poor coverage
    poor_coverage_threshold = 5.0  # Mean coverage < 5x
    low_breadth_threshold = 50.0   # Coverage breadth < 50%
    low_mapping_threshold = 30.0   # Mapping rate < 30%
    
    poor_coverage_samples = mapping_results[
        (mapping_results['Mean_Coverage'] < poor_coverage_threshold) |
        (mapping_results['Coverage_Breadth'] < low_breadth_threshold) |
        (mapping_results['Mapping_Rate'] < low_mapping_threshold)
    ]
    
    analysis_results['poor_coverage'] = {
        'count': len(poor_coverage_samples),
        'samples': poor_coverage_samples[['Sample', 'Strategy', 'Mapping_Rate', 'Mean_Coverage', 'Coverage_Breadth']].to_dict('records')
    }
    
    # 2. Identify strategies with consistently good coverage
    strategy_performance = mapping_results.groupby('Strategy').agg({
        'Mapping_Rate': ['mean', 'std'],
        'Mean_Coverage': ['mean', 'std'],
        'Coverage_Breadth': ['mean', 'std']
    })
    
    # Calculate coefficient of variation for consistency
    strategy_performance[('Mapping_Rate', 'cv')] = strategy_performance[('Mapping_Rate', 'std')] / strategy_performance[('Mapping_Rate', 'mean')]
    strategy_performance[('Mean_Coverage', 'cv')] = strategy_performance[('Mean_Coverage', 'std')] / strategy_performance[('Mean_Coverage', 'mean')]
    strategy_performance[('Coverage_Breadth', 'cv')] = strategy_performance[('Coverage_Breadth', 'std')] / strategy_performance[('Coverage_Breadth', 'mean')]
    
    analysis_results['strategy_performance'] = strategy_performance
    
    # 3. Coverage uniformity analysis
    mapping_results['Coverage_Uniformity'] = mapping_results['Median_Coverage'] / mapping_results['Mean_Coverage']
    
    # Good uniformity is close to 1 (median ≈ mean)
    uniformity_analysis = mapping_results.groupby('Strategy')['Coverage_Uniformity'].agg(['mean', 'std', 'count'])
    
    analysis_results['uniformity'] = uniformity_analysis
    
    # 4. Identify outlier samples
    outliers = {}
    for strategy in mapping_results['Strategy'].unique():
        strategy_data = mapping_results[mapping_results['Strategy'] == strategy]
        
        # Calculate z-scores for key metrics
        strategy_data['Mapping_Rate_zscore'] = np.abs((strategy_data['Mapping_Rate'] - strategy_data['Mapping_Rate'].mean()) / strategy_data['Mapping_Rate'].std())
        strategy_data['Coverage_zscore'] = np.abs((strategy_data['Mean_Coverage'] - strategy_data['Mean_Coverage'].mean()) / strategy_data['Mean_Coverage'].std())
        
        # Identify outliers (z-score > 2)
        strategy_outliers = strategy_data[
            (strategy_data['Mapping_Rate_zscore'] > 2) | 
            (strategy_data['Coverage_zscore'] > 2)
        ]
        
        if len(strategy_outliers) > 0:
            outliers[strategy] = strategy_outliers[['Sample', 'Mapping_Rate', 'Mean_Coverage', 'Mapping_Rate_zscore', 'Coverage_zscore']].to_dict('records')
    
    analysis_results['outliers'] = outliers
    
    # Save analysis results
    save_coverage_analysis_results(analysis_results, output_dir)
    
    return analysis_results

def save_coverage_analysis_results(analysis_results, output_dir):
    """Save coverage analysis results to files."""
    
    logger.info("Saving coverage analysis results")
    
    # Create analysis text report
    report_file = Path(output_dir) / "coverage_analysis_report.txt"
    
    with open(report_file, 'w') as f:
        f.write("Coverage Pattern Analysis Report\n")
        f.write("=================================\n\n")
        
        # Poor coverage samples
        f.write("1. Samples with Poor Coverage:\n")
        f.write("-" * 30 + "\n")
        poor_samples = analysis_results['poor_coverage']
        f.write(f"Total samples with issues: {poor_samples['count']}\n\n")
        
        if poor_samples['count'] > 0:
            f.write("Problematic samples:\n")
            for sample in poor_samples['samples'][:10]:  # Show first 10
                f.write(f"  {sample['Sample']} ({sample['Strategy']}): "
                       f"Mapping={sample['Mapping_Rate']:.1f}%, "
                       f"Coverage={sample['Mean_Coverage']:.1f}x, "
                       f"Breadth={sample['Coverage_Breadth']:.1f}%\n")
            
            if poor_samples['count'] > 10:
                f.write(f"  ... and {poor_samples['count'] - 10} more\n")
        
        f.write("\n")
        
        # Strategy performance
        f.write("2. Strategy Performance Summary:\n")
        f.write("-" * 32 + "\n")
        performance = analysis_results['strategy_performance']
        
        for strategy in performance.index:
            f.write(f"\n{strategy}:\n")
            f.write(f"  Mapping Rate: {performance.loc[strategy, ('Mapping_Rate', 'mean')]:.1f}% "
                   f"(±{performance.loc[strategy, ('Mapping_Rate', 'std')]:.1f})\n")
            f.write(f"  Mean Coverage: {performance.loc[strategy, ('Mean_Coverage', 'mean')]:.1f}x "
                   f"(±{performance.loc[strategy, ('Mean_Coverage', 'std')]:.1f})\n")
            f.write(f"  Coverage Breadth: {performance.loc[strategy, ('Coverage_Breadth', 'mean')]:.1f}% "
                   f"(±{performance.loc[strategy, ('Coverage_Breadth', 'std')]:.1f})\n")
        
        # Coverage uniformity
        f.write("\n3. Coverage Uniformity Analysis:\n")
        f.write("-" * 33 + "\n")
        uniformity = analysis_results['uniformity']
        
        for strategy in uniformity.index:
            uniformity_score = uniformity.loc[strategy, 'mean']
            f.write(f"{strategy}: {uniformity_score:.3f} (closer to 1.0 is better)\n")
        
        # Outliers
        f.write("\n4. Outlier Samples:\n")
        f.write("-" * 18 + "\n")
        outliers = analysis_results['outliers']
        
        if outliers:
            for strategy, strategy_outliers in outliers.items():
                f.write(f"\n{strategy} outliers ({len(strategy_outliers)}):\n")
                for outlier in strategy_outliers[:5]:  # Show first 5
                    f.write(f"  {outlier['Sample']}: "
                           f"Mapping={outlier['Mapping_Rate']:.1f}%, "
                           f"Coverage={outlier['Mean_Coverage']:.1f}x\n")
        else:
            f.write("No significant outliers detected\n")
    
    # Save detailed CSV files
    csv_dir = Path(output_dir) / "coverage_analysis_csv"
    csv_dir.mkdir(exist_ok=True)
    
    # Strategy performance CSV
    performance_df = analysis_results['strategy_performance']
    performance_df.to_csv(csv_dir / "strategy_performance.csv")
    
    # Poor coverage samples CSV
    if analysis_results['poor_coverage']['samples']:
        poor_df = pd.DataFrame(analysis_results['poor_coverage']['samples'])
        poor_df.to_csv(csv_dir / "poor_coverage_samples.csv", index=False)
    
    # Uniformity CSV
    uniformity_df = analysis_results['uniformity']
    uniformity_df.to_csv(csv_dir / "coverage_uniformity.csv")
    
    logger.info(f"Analysis report saved to {report_file}")
    logger.info(f"Detailed CSV files saved to {csv_dir}")

def generate_recommendations(mapping_results, analysis_results, output_dir):
    """Generate assembly strategy recommendations based on coverage analysis."""
    
    logger.info("Generating assembly strategy recommendations")
    
    # Calculate overall quality scores for each strategy
    strategy_scores = {}
    
    for strategy in mapping_results['Strategy'].unique():
        strategy_data = mapping_results[mapping_results['Strategy'] == strategy]
        
        # Calculate weighted quality score
        mapping_weight = 0.4
        coverage_weight = 0.3
        breadth_weight = 0.3
        
        score = (
            strategy_data['Mapping_Rate'].mean() * mapping_weight +
            np.minimum(strategy_data['Mean_Coverage'].mean(), 50) * 2 * coverage_weight +  # Cap coverage at 50x
            strategy_data['Coverage_Breadth'].mean() * breadth_weight
        )
        
        # Penalty for high variability
        cv_penalty = (
            strategy_data['Mapping_Rate'].std() / strategy_data['Mapping_Rate'].mean() +
            strategy_data['Coverage_Breadth'].std() / strategy_data['Coverage_Breadth'].mean()
        ) * 10
        
        final_score = score - cv_penalty
        
        strategy_scores[strategy] = {
            'score': final_score,
            'mapping_rate': strategy_data['Mapping_Rate'].mean(),
            'mean_coverage': strategy_data['Mean_Coverage'].mean(),
            'coverage_breadth': strategy_data['Coverage_Breadth'].mean(),
            'consistency': 100 - cv_penalty  # Higher is better
        }
    
    # Rank strategies
    ranked_strategies = sorted(strategy_scores.items(), key=lambda x: x[1]['score'], reverse=True)
    
    # Generate recommendations
    recommendations_file = Path(output_dir) / "assembly_strategy_recommendations.txt"
    
    with open(recommendations_file, 'w') as f:
        f.write("Assembly Strategy Recommendations Based on Coverage Analysis\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("Strategy Ranking (Best to Worst):\n")
        f.write("-" * 34 + "\n")
        
        for i, (strategy, metrics) in enumerate(ranked_strategies, 1):
            f.write(f"\n{i}. {strategy}\n")
            f.write(f"   Overall Score: {metrics['score']:.1f}\n")
            f.write(f"   Mapping Rate: {metrics['mapping_rate']:.1f}%\n")
            f.write(f"   Mean Coverage: {metrics['mean_coverage']:.1f}x\n")
            f.write(f"   Coverage Breadth: {metrics['coverage_breadth']:.1f}%\n")
            f.write(f"   Consistency: {metrics['consistency']:.1f}/100\n")
        
        f.write("\nRecommendations:\n")
        f.write("-" * 16 + "\n")
        
        best_strategy = ranked_strategies[0][0]
        best_score = ranked_strategies[0][1]['score']
        
        f.write(f"\n1. RECOMMENDED STRATEGY: {best_strategy}\n")
        f.write(f"   This strategy achieved the highest overall quality score ({best_score:.1f})\n")
        
        # Specific recommendations based on patterns
        poor_coverage_count = analysis_results['poor_coverage']['count']
        total_samples = len(mapping_results)
        
        if poor_coverage_count / total_samples > 0.2:
            f.write(f"\n2. COVERAGE CONCERNS:\n")
            f.write(f"   {poor_coverage_count} out of {total_samples} samples show poor coverage patterns.\n")
            f.write(f"   Consider:\n")
            f.write(f"   - Increasing sequencing depth\n")
            f.write(f"   - Checking for contamination or adapter sequences\n")
            f.write(f"   - Using more aggressive co-assembly strategies\n")
        
        # Strategy-specific recommendations
        f.write(f"\n3. STRATEGY-SPECIFIC NOTES:\n")
        
        for strategy, metrics in strategy_scores.items():
            if metrics['mapping_rate'] < 40:
                f.write(f"   - {strategy}: Low mapping rate ({metrics['mapping_rate']:.1f}%) suggests assembly fragmentation\n")
            
            if metrics['mean_coverage'] < 5:
                f.write(f"   - {strategy}: Low coverage ({metrics['mean_coverage']:.1f}x) may affect assembly quality\n")
            
            if metrics['coverage_breadth'] < 60:
                f.write(f"   - {strategy}: Low breadth ({metrics['coverage_breadth']:.1f}%) indicates incomplete assemblies\n")
            
            if metrics['consistency'] < 70:
                f.write(f"   - {strategy}: High variability (consistency: {metrics['consistency']:.1f}) suggests inconsistent performance\n")
        
        f.write(f"\n4. NEXT STEPS:\n")
        f.write(f"   - Use {best_strategy} for downstream analyses\n")
        f.write(f"   - Validate results with CheckV quality assessment\n")
        f.write(f"   - Consider combining strategies if appropriate\n")
        f.write(f"   - Investigate samples flagged for poor coverage\n")
    
    logger.info(f"Recommendations saved to {recommendations_file}")
    
    return ranked_strategies

def main():
    parser = argparse.ArgumentParser(description="Analyze coverage patterns from read mapping results")
    parser.add_argument("--mapping-dir", default="results/quality_assessment/read_mapping",
                       help="Directory containing read mapping results")
    parser.add_argument("--output-dir", default="results/quality_assessment/coverage_analysis",
                       help="Output directory for coverage analysis")
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Load mapping results
        mapping_results, individual_results = load_mapping_results(args.mapping_dir)
        
        # Calculate coverage statistics
        coverage_stats = calculate_coverage_distributions(mapping_results)
        
        # Save coverage statistics
        coverage_stats.to_csv(output_dir / "coverage_statistics.csv")
        logger.info(f"Coverage statistics saved to {output_dir / 'coverage_statistics.csv'}")
        
        # Create visualizations
        create_coverage_visualizations(mapping_results, output_dir)
        
        # Analyze coverage patterns
        analysis_results = analyze_coverage_patterns(mapping_results, output_dir)
        
        # Generate recommendations
        recommendations = generate_recommendations(mapping_results, analysis_results, output_dir)
        
        logger.info("Coverage analysis completed successfully!")
        logger.info(f"Results saved to: {output_dir}")
        
        # Print summary
        print("\nCoverage Analysis Summary:")
        print("=" * 26)
        print(f"Total samples analyzed: {len(mapping_results)}")
        print(f"Strategies compared: {len(mapping_results['Strategy'].unique())}")
        print(f"Samples with poor coverage: {analysis_results['poor_coverage']['count']}")
        print(f"\nTop strategy: {recommendations[0][0]} (score: {recommendations[0][1]['score']:.1f})")
        print(f"Results saved to: {output_dir}")
        
    except Exception as e:
        logger.error(f"Coverage analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()