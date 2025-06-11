# Troubleshooting Guide

This guide provides solutions to common issues encountered when using the Viral Metagenomic Assembly Toolkit.

## Installation Issues

### Missing Dependencies

**Problem**: Script fails with "command not found" errors

**Solution**:
```bash
# Check which tools are missing
./scripts/run_pipeline.sh --help

# Install missing tools using conda
conda install -c bioconda sourmash megahit spades checkv bwa samtools bedtools seqtk

# Or using specific package managers
# For Sourmash
pip install sourmash

# For CheckV
checkv download_database /path/to/checkv-db
```

### Python/R Package Issues

**Problem**: Import errors for Python or R packages

**Solution**:
```bash
# Python packages
pip install pandas numpy matplotlib seaborn scikit-learn biopython

# R packages
Rscript -e "install.packages(c('vegan', 'ggplot2', 'dplyr', 'viridis', 'RColorBrewer'))"
```

## Data Preparation Issues

### Input File Format Problems

**Problem**: "No samples found" or file format errors

**Solution**:
1. Check file naming convention matches config:
   ```
   sample_R1.fastq.gz
   sample_R2.fastq.gz
   ```

2. Verify file paths in config.yaml:
   ```yaml
   input:
     reads_directory: "data/cleaned_reads"
     read_pattern:
       r1: "{sample}_R1.fastq.gz"
       r2: "{sample}_R2.fastq.gz"
   ```

3. Check file permissions:
   ```bash
   ls -la data/cleaned_reads/
   chmod 644 data/cleaned_reads/*.fastq.gz
   ```

### Metadata File Issues

**Problem**: PERMANOVA analysis fails due to metadata issues

**Solution**:
1. Check CSV format (comma-separated, proper headers)
2. Ensure Sample column matches FASTQ file names
3. Remove special characters from column names
4. Check for missing values:
   ```bash
   head -5 examples/metadata_template.csv
   ```

## Assembly Issues

### Memory Errors

**Problem**: Assembly fails with out-of-memory errors

**Solutions**:
1. Reduce memory usage:
   ```yaml
   resources:
     memory_gb: 16  # Reduce from 32
   
   assembly:
     assemblers: ["megahit"]  # Use only MEGAHIT (more memory efficient)
   ```

2. Reduce co-assembly group sizes:
   ```yaml
   kmer_analysis:
     max_group_size: 4  # Reduce from 8
   ```

3. Use individual assemblies only:
   ```yaml
   assembly:
     strategies:
       individual: true
       strategic_coassembly: false
       global_coassembly: false
   ```

### Assembly Quality Issues

**Problem**: Very low N50 or poor assembly statistics

**Solutions**:
1. Check input data quality:
   ```bash
   fastqc data/cleaned_reads/*.fastq.gz
   ```

2. Adjust assembly parameters:
   ```yaml
   assembly:
     min_contig_length: 300  # Reduce from 500
     megahit:
       presets: "meta-large"  # Try different preset
   ```

3. Check for contamination:
   ```bash
   # Verify host filtering was effective
   kraken2 --db standard data/cleaned_reads/sample_R1.fastq.gz
   ```

### Disk Space Issues

**Problem**: Pipeline fails due to insufficient disk space

**Solutions**:
1. Clean up intermediate files:
   ```bash
   rm -rf results/assemblies/*/temp/
   rm -rf results/*/megahit_*/intermediate_contigs/
   ```

2. Use compression:
   ```yaml
   output:
     compress_outputs: true
     keep_intermediate: false
   ```

3. Redirect temp files:
   ```yaml
   resources:
     temp_dir: "/scratch/tmp"  # Use larger temp space
   ```

## Quality Assessment Issues

### CheckV Database Problems

**Problem**: CheckV fails with database errors

**Solutions**:
1. Download/update CheckV database:
   ```bash
   checkv download_database /opt/checkv-db-v1.5
   ```

2. Update config with correct path:
   ```yaml
   databases:
     checkv: "/opt/checkv-db-v1.5"
   ```

3. Check database permissions:
   ```bash
   ls -la /opt/checkv-db-v1.5/
   ```

### Low Mapping Rates

**Problem**: Very low read mapping rates (<30%)

**Solutions**:
1. Check for adapter contamination:
   ```bash
   fastqc data/cleaned_reads/*.fastq.gz
   # Look for adapter content
   ```

2. Verify host filtering:
   ```bash
   # Check if host sequences remain
   bowtie2 -x human_genome -1 sample_R1.fastq -2 sample_R2.fastq --un-conc unmapped.fastq
   ```

3. Adjust mapping parameters:
   ```yaml
   quality_assessment:
     mapping:
       min_mapq: 5  # Reduce from 10
   ```

### No High-Quality Viral Contigs

**Problem**: CheckV finds no high-quality viral contigs

**Solutions**:
1. Check contig length distribution:
   ```bash
   seqkit stats results/assemblies/*/contigs.fasta
   ```

2. Lower minimum length threshold:
   ```yaml
   quality_assessment:
     min_contig_length_checkv: 500  # Reduce from 1000
   ```

3. Verify viral enrichment:
   ```bash
   # Check for viral signatures
   hmmsearch --tblout viral_genes.txt /path/to/viral_hmms contigs.fasta
   ```

## Performance Issues

### Slow K-mer Analysis

**Problem**: Sourmash analysis takes too long

**Solutions**:
1. Reduce sketch size:
   ```yaml
   kmer_analysis:
     sketch_size: 5000  # Reduce from 10000
   ```

2. Use larger k-mer size:
   ```yaml
   kmer_analysis:
     kmer_size: 51  # Increase from 31
   ```

### Pipeline Hangs

**Problem**: Pipeline stops responding

**Solutions**:
1. Check system resources:
   ```bash
   top
   df -h
   ```

2. Reduce parallel jobs:
   ```yaml
   resources:
     max_parallel_jobs: 2  # Reduce from 4
   ```

3. Resume from checkpoint:
   ```bash
   RESUME=true ./scripts/run_pipeline.sh
   ```

## Analysis Issues

### PERMANOVA Fails

**Problem**: Variable analysis produces errors

**Solutions**:
1. Check sample size:
   - Need at least 10 samples for reliable PERMANOVA
   - Ensure sufficient replication within groups

2. Simplify metadata:
   ```R
   # Remove columns with too many unique values
   # Keep only essential variables
   ```

3. Handle missing values:
   ```R
   # Remove samples with missing critical metadata
   complete_cases <- complete.cases(metadata[,essential_vars])
   ```

### No Significant Variables

**Problem**: PERMANOVA finds no significant variables

**Solutions**:
1. This is often expected - not all datasets have strong metadata effects
2. Use k-mer similarity as primary grouping criterion
3. Consider additional metadata variables
4. Check for batch effects

### Comparison Results Look Wrong

**Problem**: Strategy rankings seem counterintuitive

**Solutions**:
1. Check input data quality across all strategies
2. Verify all assembly strategies completed successfully
3. Review individual quality metrics:
   ```bash
   cat results/comparison/integrated_comparison_results.csv
   ```

4. Consider dataset-specific factors

## Hardware-Specific Issues

### macOS Issues

**Problem**: Tools behave differently on macOS

**Solutions**:
1. Use GNU versions of tools:
   ```bash
   brew install coreutils gnu-sed
   export PATH="/usr/local/opt/gnu-sed/libexec/gnubin:$PATH"
   ```

2. Install XCode command line tools:
   ```bash
   xcode-select --install
   ```

### Windows/WSL Issues

**Problem**: Path or permission issues on Windows

**Solutions**:
1. Use WSL2 for better compatibility
2. Ensure UNIX line endings:
   ```bash
   dos2unix scripts/*.sh
   ```

3. Fix permissions:
   ```bash
   chmod +x scripts/*.sh scripts/*/*.sh
   ```

## Advanced Troubleshooting

### Debug Mode

Enable detailed logging:
```bash
export LOG_LEVEL=DEBUG
./scripts/run_pipeline.sh config.yaml 2>&1 | tee pipeline.log
```

### Manual Step Execution

Run individual steps for debugging:
```bash
# Test k-mer analysis only
bash scripts/kmer_profiling/01_generate_sketches.sh

# Test specific assembly
bash scripts/assembly/individual_assembly.sh
```

### Memory Profiling

Monitor memory usage:
```bash
# Track memory during assembly
/usr/bin/time -v bash scripts/assembly/individual_assembly.sh
```

### Validate Dependencies

Check all tool versions:
```bash
sourmash --version
megahit --version
checkv --version
python --version
R --version
```

## Getting Help

If issues persist:

1. **Check logs**: Look in `results/*/logs/` for detailed error messages
2. **Simplify dataset**: Test with a smaller subset of samples
3. **Update software**: Ensure all tools are current versions
4. **Community support**: Post issues on GitHub with:
   - Error messages
   - Configuration file
   - System information
   - Steps to reproduce

## Common Error Messages

### "No module named 'pandas'"
**Solution**: Install Python dependencies
```bash
pip install pandas numpy matplotlib seaborn
```

### "Error in library(vegan)"
**Solution**: Install R packages
```bash
Rscript -e "install.packages('vegan')"
```

### "Permission denied"
**Solution**: Fix file permissions
```bash
chmod +x scripts/*.sh
chmod 644 data/*.fastq.gz
```

### "Disk quota exceeded"
**Solution**: Clean up space or use different temp directory
```bash
df -h  # Check disk usage
rm -rf /tmp/megahit_*  # Clean temp files
```

### "Cannot allocate memory"
**Solution**: Reduce memory requirements or use smaller datasets
```yaml
resources:
  memory_gb: 16
assembly:
  assemblers: ["megahit"]  # More memory efficient
```