# Viral Metagenomic Assembly Optimization Toolkit

A comprehensive toolkit for optimizing viral metagenomic assembly strategies, particularly for virus-like particle (VLP) sequencing data. This repository provides standardized protocols and scripts to determine the optimal assembly approach for your viral metagenomic datasets.

## Table of Contents

- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Quick Start](#quick-start)
- [Standard Operating Procedure](#standard-operating-procedure)
- [Directory Structure](#directory-structure)
- [Computational Requirements](#computational-requirements)
- [Troubleshooting](#troubleshooting)
- [Integration with Downstream Analyses](#integration-with-downstream-analyses)

## Introduction

Viral metagenomic assembly presents unique challenges due to:
- High genetic diversity within viral populations
- Uneven coverage distributions
- Strain-level variations
- Small genome sizes and fragmented assemblies
- Contamination from host DNA and other microorganisms

This toolkit addresses these challenges by providing a systematic approach to:
1. **Profile sample similarity** using k-mer analysis
2. **Identify important grouping variables** (patient ID, geography, time points)
3. **Execute multiple assembly strategies** (individual, strategic co-assembly, global co-assembly, meta-assembly)
4. **Assess assembly quality** comprehensively
5. **Compare and select optimal approaches** based on quantitative metrics

## Prerequisites

### Essential Data Preprocessing Steps

Before using this toolkit, ensure your sequencing data has undergone:

1. **Adapter Removal**: Remove sequencing adapters and low-quality bases

2. **Host Filtering**: Remove host contamination 

3. **PCR Bias Correction**: Remove PCR duplicates
   ```bash
   # Example using BBTools
   clumpify.sh in1=sample_R1.fastq in2=sample_R2.fastq out1=sample_R1_dedup.fastq out2=sample_R2_dedup.fastq dedupe
   ```

### Required Software

- **Sourmash** (v4.0+): K-mer profiling and comparison
- **MEGAHIT** (v1.2+): Fast metagenomic assembler
- **metaSPAdes** (v3.15+): Metagenomic assembler
- **CheckV** (v1.0+): Viral contig quality assessment
- **BWA** or **Bowtie2**: Read mapping
- **R** (v4.0+) with packages: vegan, ggplot2, dplyr, phyloseq
- **Python** (v3.8+) with packages: pandas, numpy, matplotlib, seaborn, scikit-learn

## Quick Start

```bash
# Clone the repository
git clone https://github.com/megjohnson1999/virome-assembly-optimizer.git
cd virome-assembly-optimizer

# Set up configuration
cp config/config_template.yaml config/config.yaml
# Edit config.yaml with your specific parameters

# Run the complete pipeline
bash scripts/run_pipeline.sh config/config.yaml
```

## Standard Operating Procedure

### Step 1: K-mer Profiling for Sample Similarity Assessment

**Objective**: Identify samples with similar viral compositions that may benefit from co-assembly.

```bash
# Generate k-mer sketches for all samples
bash scripts/kmer_profiling/01_generate_sketches.sh

# Compare samples and create distance matrix
python scripts/kmer_profiling/02_compare_samples.py

# Visualize sample relationships
Rscript scripts/kmer_profiling/03_visualize_similarity.R
```

**Key Decision Points**:
- Identify samples with high k-mer similarity as candidates for co-assembly
- Consider batch effects and technical replicates
- Document clustering patterns for downstream grouping

### Step 2: Determine Important Grouping Variables

**Objective**: Identify which metadata variables (patient ID, geography, time point) significantly influence viral community structure.

```bash
# Prepare metadata file (see examples/metadata_template.csv)
# Run PERMANOVA analysis
Rscript scripts/variable_analysis/01_permanova_analysis.R

# Visualize variable importance
Rscript scripts/variable_analysis/02_plot_variable_importance.R
```

**Interpretation Guidelines**:
- Variables with significant p-values and meaningful effect sizes are considered important
- Prioritize biological variables over technical ones
- Use results to design co-assembly groups

### Step 3: Execute Multiple Assembly Strategies

All strategies produce a final community-level assembly for fair comparison.

**Strategy A: Individual + Meta-assembly**
```bash
# Step 1: Individual sample assemblies
bash scripts/assembly/individual_assembly.sh
# Step 2: Meta-assembly to create community assembly
bash scripts/assembly/meta_assembly.sh --input-type individual
```

**Strategy B: Strategic Co-assemblies + Meta-assembly** (based on similarity and important variables)
```bash
# Step 1: Create co-assembly groups based on analysis results
python scripts/assembly/create_coassembly_groups.py
# Step 2: Run strategic co-assemblies
bash scripts/assembly/strategic_coassembly.sh
# Step 3: Meta-assembly to create community assembly
bash scripts/assembly/meta_assembly.sh --input-type strategic
```

**Strategy C: Global Co-assembly** (all samples together - already community-level)
```bash
bash scripts/assembly/global_coassembly.sh
```

### Step 4: Quality Assessment

**CheckV Analysis**:
```bash
bash scripts/quality_assessment/01_run_checkv.sh
```

**Read Mapping and Coverage Analysis**:
```bash
bash scripts/quality_assessment/02_read_mapping.sh
python scripts/quality_assessment/03_analyze_coverage.py
```

**Contig Statistics**:
```bash
python scripts/quality_assessment/04_contig_stats.py
```

### Step 5: Strategy Comparison and Selection

```bash
# Compare assembly strategies
python scripts/comparison/01_compare_strategies.py

# Generate comprehensive report
Rscript scripts/comparison/02_generate_report.R
```

**Selection Criteria**:
1. **Viral contig count**: Higher number of high-quality viral contigs
2. **Completeness**: Proportion of complete viral genomes
3. **Contamination**: Lower host/bacterial contamination
4. **Coverage uniformity**: More even read coverage
5. **Gene recovery**: Recovery of essential viral genes

**Important Note on Thresholds**: This toolkit provides data distributions and statistical summaries to guide threshold selection. Users should define specific thresholds (e.g., similarity cutoffs, coverage requirements) based on their dataset characteristics, sequencing depth, and research objectives rather than using arbitrary predetermined values.

## Directory Structure

```
virome-assembly-toolkit/
├── README.md
├── config/
│   ├── config_template.yaml
│   └── checkv_config.yaml
├── scripts/
│   ├── kmer_profiling/
│   │   ├── 01_generate_sketches.sh
│   │   ├── 02_compare_samples.py
│   │   └── 03_visualize_similarity.R
│   ├── variable_analysis/
│   │   ├── 01_permanova_analysis.R
│   │   └── 02_plot_variable_importance.R
│   ├── assembly/
│   │   ├── individual_assembly.sh
│   │   ├── strategic_coassembly.sh
│   │   ├── global_coassembly.sh
│   │   ├── meta_assembly.sh
│   │   └── create_coassembly_groups.py
│   ├── quality_assessment/
│   │   ├── 01_run_checkv.sh
│   │   ├── 02_read_mapping.sh
│   │   ├── 03_analyze_coverage.py
│   │   └── 04_contig_stats.py
│   ├── comparison/
│   │   ├── 01_compare_strategies.py
│   │   └── 02_generate_report.R
│   └── run_pipeline.sh
├── examples/
│   ├── metadata_template.csv
│   ├── sample_list.txt
│   └── expected_output/
├── data/
│   └── (user data directory)
└── docs/
    ├── troubleshooting.md
    └── computational_requirements.md
```

## Computational Requirements

### Minimum Requirements
- **CPU**: 8 cores
- **RAM**: 32 GB
- **Storage**: 100 GB free space per dataset
- **Time**: 2-24 hours depending on dataset size

### Recommended Resources
- **CPU**: 16-32 cores
- **RAM**: 64-128 GB
- **Storage**: 500 GB SSD
- **Time**: 4-12 hours for most datasets

### Scaling Guidelines
- **Individual assemblies**: RAM requirements scale with sample size
- **Co-assemblies**: RAM scales with combined sample size
- **Global co-assembly**: May require substantial RAM for large datasets

## Troubleshooting

### Common Issues

**Low k-mer similarity between samples**:
- Check for batch effects or technical issues
- Consider individual assemblies or smaller co-assembly groups
- Verify preprocessing quality

**Assembly failures**:
- Reduce k-mer size for low-coverage samples
- Adjust memory settings for large co-assemblies
- Check input file integrity

**Poor CheckV results**:
- Verify CheckV database is current
- Check for host contamination in assemblies
- Consider adjusting contig length thresholds

**High memory usage**:
- Process samples in smaller batches
- Use memory-efficient assemblers (MEGAHIT vs metaSPAdes)
- Implement temporary file cleanup

For detailed troubleshooting, see [docs/troubleshooting.md](docs/troubleshooting.md).

## Integration with Downstream Analyses

### Viral Taxonomy Assignment
```bash
# Example using vConTACT2
vcontact2 --raw-proteins viral_proteins.faa --rel-mode 'Diamond' --proteins-fp gene_to_genome.csv --c1-bin /path/to/cluster_one --output-dir vcontact2_output
```

### Functional Annotation
```bash
# Example using DRAM-v
DRAM-v.py annotate -i viral_contigs.fasta -o dram_annotation --threads 16
```

### Phylogenetic Analysis
```bash
# Example using IQ-TREE
iqtree -s concatenated_alignment.fasta -bb 1000 -nt AUTO
```

### Host Prediction
```bash
# Example using iPHoP
iphop predict --fa_file viral_contigs.fasta --db_dir /path/to/iphop_db --out_dir host_prediction
```

## Contributing

We welcome contributions! Please see our [contributing guidelines](CONTRIBUTING.md) for details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

For questions and support:
- Open an issue on GitHub
