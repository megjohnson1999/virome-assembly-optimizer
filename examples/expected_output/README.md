# Expected Output Examples

This directory contains example outputs from the viral metagenomic assembly toolkit to help users understand what to expect from each analysis step.

## Directory Structure

When you run the complete pipeline, you can expect the following output structure:

```
results/
├── kmer_sketches/                    # Sourmash k-mer sketches
├── similarity_analysis/              # K-mer similarity results
│   ├── similarity_matrix.csv
│   ├── pairwise_similarities.csv
│   ├── clustering_suggestions.json
│   └── plots/
├── variable_analysis/                # PERMANOVA results
│   ├── permanova_results.csv
│   ├── significant_variables.txt
│   ├── assembly_recommendations.txt
│   └── plots/
├── coassembly_groups/               # Strategic grouping
│   ├── coassembly_groups.json
│   ├── assembly_summary.txt
│   └── commands/
├── assemblies/                      # Assembly outputs
│   ├── individual/
│   ├── strategic_coassembly/
│   ├── global_coassembly/
│   └── meta_assembly/
├── quality_assessment/              # Quality analysis
│   ├── checkv/
│   ├── read_mapping/
│   ├── coverage_analysis/
│   └── contig_stats/
└── comparison/                      # Final comparison
    ├── integrated_comparison_results.csv
    ├── final_assembly_recommendations.txt
    └── publication_figures/
```

## Key Output Files

### 1. Assembly Strategy Recommendations
- **File**: `results/comparison/final_assembly_recommendations.txt`
- **Description**: Comprehensive recommendations for optimal assembly strategy
- **Use**: Primary decision-making document

### 2. Quality Comparison Matrix
- **File**: `results/comparison/integrated_comparison_results.csv`
- **Description**: Detailed metrics for all assembly strategies
- **Use**: Quantitative comparison of all approaches

### 3. K-mer Similarity Results
- **File**: `results/similarity_analysis/clustering_suggestions.json`
- **Description**: Suggested sample groupings based on k-mer similarity
- **Use**: Guides strategic co-assembly decisions

### 4. Variable Importance Results
- **File**: `results/variable_analysis/permanova_results.csv`
- **Description**: Statistical significance of metadata variables
- **Use**: Identifies which factors influence viral community structure

### 5. CheckV Quality Assessment
- **File**: `results/quality_assessment/checkv/combined/checkv_comparison.csv`
- **Description**: Viral genome completeness and quality metrics
- **Use**: Evaluates viral genome recovery success

### 6. Publication-Ready Figures
- **Directory**: `results/comparison/publication_figures/`
- **Description**: High-resolution figures for publication
- **Use**: Direct inclusion in manuscripts and presentations

## Expected Runtime

Typical runtimes for different dataset sizes:

| Dataset Size | Samples | Runtime | Memory | Storage |
|--------------|---------|---------|---------|---------|
| Small        | 5-10    | 2-4 hours | 16 GB | 50 GB |
| Medium       | 15-30   | 6-12 hours | 32 GB | 200 GB |
| Large        | 50+     | 1-2 days | 64+ GB | 500+ GB |

## Quality Thresholds

Expected quality ranges for good assemblies:

- **N50**: > 5,000 bp (excellent: > 10,000 bp)
- **Mapping Rate**: > 70% (excellent: > 85%)
- **High-Quality Viral Contigs**: > 20% (excellent: > 40%)
- **Coverage Breadth**: > 60% (excellent: > 80%)
- **Contamination Rate**: < 10% (excellent: < 5%)

## Troubleshooting Common Issues

### Low Assembly Quality
- Check input data quality (adapter removal, host filtering)
- Consider different assembly parameters
- Verify sufficient sequencing depth

### High Memory Usage
- Reduce co-assembly group sizes
- Use MEGAHIT instead of metaSPAdes
- Process samples in smaller batches

### Poor Mapping Rates
- Check for adapter contamination
- Verify host filtering effectiveness
- Consider assembly fragmentation issues

### No High-Quality Viral Contigs
- Check CheckV database version
- Verify input data is viral-enriched
- Consider longer contigs (increase min length thresholds)

## Validation Recommendations

1. **Cross-validation**: Compare results with alternative tools
2. **Manual inspection**: Review top contigs individually
3. **Functional annotation**: Confirm viral genes in high-quality contigs
4. **Phylogenetic analysis**: Validate taxonomic assignments
5. **Experimental validation**: PCR/qPCR confirmation when possible