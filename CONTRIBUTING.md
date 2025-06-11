# Contributing to Viral Metagenomic Assembly Toolkit

Thank you for your interest in contributing to the Viral Metagenomic Assembly Toolkit! This project aims to provide a comprehensive, scientifically rigorous approach to optimizing viral metagenomic assembly strategies.

## How to Contribute

### Reporting Issues

If you encounter bugs, have feature requests, or find areas for improvement:

1. **Search existing issues** to avoid duplicates
2. **Open a new issue** with:
   - Clear, descriptive title
   - Detailed description of the problem/request
   - Steps to reproduce (for bugs)
   - System information (OS, software versions)
   - Example data/config files (if applicable)

### Contributing Code

1. **Fork the repository**
2. **Create a feature branch**: `git checkout -b feature/your-feature-name`
3. **Make your changes** following our coding standards
4. **Test thoroughly** with example datasets
5. **Update documentation** as needed
6. **Submit a pull request**

## Development Guidelines

### Code Standards

**Shell Scripts:**
- Use `#!/bin/bash` with `set -euo pipefail`
- Include clear comments and function documentation
- Follow consistent variable naming (UPPER_CASE for globals)
- Add error handling and validation

**Python Scripts:**
- Follow PEP 8 style guidelines
- Include docstrings for all functions
- Use type hints where appropriate
- Add comprehensive error handling
- Include logging with appropriate levels

**R Scripts:**
- Follow tidyverse style guide
- Document all functions with roxygen2-style comments
- Use consistent naming conventions
- Include package dependency checks

### Testing

Before submitting changes:

1. **Test with example data**:
   ```bash
   # Test complete pipeline
   ./scripts/run_pipeline.sh config/config_template.yaml
   ```

2. **Validate script execution**:
   ```bash
   # Test individual components
   bash scripts/kmer_profiling/01_generate_sketches.sh
   python scripts/quality_assessment/03_analyze_coverage.py --help
   ```

3. **Check cross-platform compatibility** (Linux/macOS)

### Documentation

When contributing:
- Update README.md if adding new features
- Add/update docstrings and comments
- Include example usage in commit messages
- Update troubleshooting guide for new issues

## Specific Contribution Areas

### High Priority
- **New assembly strategies** (e.g., hybrid approaches)
- **Additional quality metrics** (e.g., gene completeness)
- **Performance optimizations** (memory, speed)
- **Cloud computing integration** (AWS, Google Cloud)

### Medium Priority  
- **Additional assemblers** (Flye, Unicycler, etc.)
- **Containerization** (Docker, Singularity)
- **Workflow managers** (Snakemake, Nextflow)
- **Database integrations** (NCBI, ENA)

### Documentation & Usability
- **Tutorial development** (step-by-step guides)
- **Video tutorials** (screencasts)
- **Benchmark datasets** (standard test cases)
- **Performance benchmarking** (systematic comparisons)

## Scientific Contributions

### Algorithm Improvements
- **Novel similarity metrics** for sample grouping
- **Advanced statistical methods** for variable analysis
- **Machine learning approaches** for strategy selection
- **Quality prediction models**

### Validation Studies
- **Systematic benchmarking** across diverse datasets
- **Comparison with other tools** (comprehensive evaluation)
- **Parameter optimization** studies
- **Cross-validation** approaches

## Code Review Process

All contributions go through peer review:

1. **Automated checks** (syntax, style, basic functionality)
2. **Scientific review** (algorithmic correctness, statistical validity)
3. **Usability testing** (documentation clarity, ease of use)
4. **Performance evaluation** (computational efficiency)

## Community Guidelines

### Be Respectful
- Welcome newcomers and different perspectives
- Provide constructive feedback
- Acknowledge contributions from others
- Follow professional communication standards

### Be Collaborative
- Share knowledge and expertise
- Help others learn and improve
- Coordinate on overlapping efforts
- Document decisions and rationale

### Be Scientific
- Prioritize accuracy and reproducibility
- Cite relevant literature and methods
- Validate claims with appropriate testing
- Consider statistical significance and effect sizes

## Getting Started

### Development Environment Setup

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/virome-assembly-toolkit.git
   cd virome-assembly-toolkit
   ```

2. **Install dependencies**:
   ```bash
   # Using conda (recommended)
   conda env create -f environment.yml
   conda activate virome-toolkit
   
   # Or install manually
   conda install -c bioconda sourmash megahit spades checkv bwa samtools
   pip install pandas numpy matplotlib seaborn biopython
   ```

3. **Test installation**:
   ```bash
   ./scripts/run_pipeline.sh --help
   ```

### First Contribution Ideas

**Good first issues:**
- Fix typos in documentation
- Add error messages for common mistakes  
- Improve log output formatting
- Add parameter validation
- Create additional test cases

**Medium complexity:**
- Add new visualization options
- Implement additional file format support
- Optimize existing algorithms
- Add progress bars/status indicators

## Release Process

### Version Numbering
We follow semantic versioning (MAJOR.MINOR.PATCH):
- **MAJOR**: Incompatible API changes
- **MINOR**: New functionality (backward compatible)
- **PATCH**: Bug fixes (backward compatible)

### Release Checklist
- [ ] All tests pass
- [ ] Documentation updated
- [ ] CHANGELOG.md updated
- [ ] Version numbers updated
- [ ] Release notes prepared

## Questions?

- **General questions**: Open a GitHub issue with the "question" label
- **Development discussion**: Use GitHub Discussions
- **Private inquiries**: Contact maintainers directly

## Acknowledgments

Contributors will be acknowledged in:
- README.md contributor section
- Publication acknowledgments (for significant contributions)
- Release notes
- Documentation credits

Thank you for helping make viral metagenomic analysis more accessible and reliable!