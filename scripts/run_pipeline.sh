#!/bin/bash

# Master pipeline script for viral metagenomic assembly optimization
# Executes the complete workflow from k-mer analysis to final recommendations

set -euo pipefail

# Default configuration
CONFIG_FILE=${1:-"config/config.yaml"}
RESUME=${RESUME:-false}
FORCE=${FORCE:-false}

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    local color=$1
    local message=$2
    echo -e "${color}[$(date '+%Y-%m-%d %H:%M:%S')] ${message}${NC}"
}

print_info() { print_status "$BLUE" "INFO: $1"; }
print_success() { print_status "$GREEN" "SUCCESS: $1"; }
print_warning() { print_status "$YELLOW" "WARNING: $1"; }
print_error() { print_status "$RED" "ERROR: $1"; }

# Function to check if a step is completed
check_step_completion() {
    local step_name=$1
    local marker_file="results/.${step_name}_completed"
    
    if [[ -f "$marker_file" && "$FORCE" != "true" ]]; then
        return 0  # Step completed
    else
        return 1  # Step not completed
    fi
}

# Function to mark step as completed
mark_step_completed() {
    local step_name=$1
    local marker_file="results/.${step_name}_completed"
    
    mkdir -p results
    touch "$marker_file"
    echo "$(date)" > "$marker_file"
}

# Function to validate configuration
validate_config() {
    print_info "Validating configuration file: $CONFIG_FILE"
    
    if [[ ! -f "$CONFIG_FILE" ]]; then
        print_error "Configuration file not found: $CONFIG_FILE"
        print_info "Copy config/config_template.yaml to $CONFIG_FILE and customize it"
        exit 1
    fi
    
    # Check for required tools
    local required_tools=("sourmash" "python3" "Rscript" "megahit" "checkv" "bwa" "samtools")
    local missing_tools=()
    
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        print_error "Missing required tools: ${missing_tools[*]}"
        print_info "Please install missing tools before running the pipeline"
        exit 1
    fi
    
    print_success "Configuration and dependencies validated"
}

# Function to run k-mer profiling
run_kmer_profiling() {
    local step_name="kmer_profiling"
    
    if check_step_completion "$step_name"; then
        print_info "K-mer profiling already completed, skipping..."
        return 0
    fi
    
    print_info "Step 1: Running k-mer profiling analysis"
    
    # Generate k-mer sketches
    if ! bash scripts/kmer_profiling/01_generate_sketches.sh; then
        print_error "Failed to generate k-mer sketches"
        return 1
    fi
    
    # Compare samples
    if ! python scripts/kmer_profiling/02_compare_samples.py; then
        print_error "Failed to compare k-mer samples"
        return 1
    fi
    
    # Visualize similarity
    if ! Rscript scripts/kmer_profiling/03_visualize_similarity.R; then
        print_error "Failed to create similarity visualizations"
        return 1
    fi
    
    mark_step_completed "$step_name"
    print_success "K-mer profiling completed"
}

# Function to run variable analysis
run_variable_analysis() {
    local step_name="variable_analysis"
    
    if check_step_completion "$step_name"; then
        print_info "Variable analysis already completed, skipping..."
        return 0
    fi
    
    print_info "Step 2: Running variable importance analysis"
    
    # Run PERMANOVA analysis
    if ! Rscript scripts/variable_analysis/01_permanova_analysis.R; then
        print_error "Failed to run PERMANOVA analysis"
        return 1
    fi
    
    # Create importance visualizations
    if ! Rscript scripts/variable_analysis/02_plot_variable_importance.R; then
        print_error "Failed to create variable importance plots"
        return 1
    fi
    
    mark_step_completed "$step_name"
    print_success "Variable analysis completed"
}

# Function to create co-assembly groups
create_coassembly_groups() {
    local step_name="coassembly_groups"
    
    if check_step_completion "$step_name"; then
        print_info "Co-assembly groups already created, skipping..."
        return 0
    fi
    
    print_info "Step 3: Creating strategic co-assembly groups"
    
    if ! python scripts/assembly/create_coassembly_groups.py; then
        print_error "Failed to create co-assembly groups"
        return 1
    fi
    
    mark_step_completed "$step_name"
    print_success "Co-assembly groups created"
}

# Function to run assemblies
run_assemblies() {
    local step_name="assemblies"
    
    if check_step_completion "$step_name"; then
        print_info "Assemblies already completed, skipping..."
        return 0
    fi
    
    print_info "Step 4: Running assembly strategies"
    
    # Individual assemblies
    print_info "  Running individual assemblies..."
    if ! bash scripts/assembly/individual_assembly.sh; then
        print_warning "Individual assemblies failed"
    fi
    
    # Strategic co-assemblies
    print_info "  Running strategic co-assemblies..."
    if ! bash scripts/assembly/strategic_coassembly.sh; then
        print_warning "Strategic co-assemblies failed"
    fi
    
    # Global co-assembly
    print_info "  Running global co-assembly..."
    if ! bash scripts/assembly/global_coassembly.sh; then
        print_warning "Global co-assembly failed"
    fi
    
    # Meta-assembly
    print_info "  Running meta-assembly..."
    if ! bash scripts/assembly/meta_assembly.sh; then
        print_warning "Meta-assembly failed"
    fi
    
    mark_step_completed "$step_name"
    print_success "Assembly strategies completed"
}

# Function to run quality assessment
run_quality_assessment() {
    local step_name="quality_assessment"
    
    if check_step_completion "$step_name"; then
        print_info "Quality assessment already completed, skipping..."
        return 0
    fi
    
    print_info "Step 5: Running quality assessment"
    
    # CheckV analysis
    print_info "  Running CheckV analysis..."
    if ! bash scripts/quality_assessment/01_run_checkv.sh; then
        print_warning "CheckV analysis failed"
    fi
    
    # Read mapping analysis
    print_info "  Running read mapping analysis..."
    if ! bash scripts/quality_assessment/02_read_mapping.sh; then
        print_warning "Read mapping analysis failed"
    fi
    
    # Coverage analysis
    print_info "  Running coverage analysis..."
    if ! python scripts/quality_assessment/03_analyze_coverage.py; then
        print_warning "Coverage analysis failed"
    fi
    
    # Contig statistics
    print_info "  Running contig statistics analysis..."
    if ! python scripts/quality_assessment/04_contig_stats.py; then
        print_warning "Contig statistics analysis failed"
    fi
    
    mark_step_completed "$step_name"
    print_success "Quality assessment completed"
}

# Function to run strategy comparison
run_strategy_comparison() {
    local step_name="strategy_comparison"
    
    if check_step_completion "$step_name"; then
        print_info "Strategy comparison already completed, skipping..."
        return 0
    fi
    
    print_info "Step 6: Running strategy comparison and generating reports"
    
    # Compare strategies
    print_info "  Comparing assembly strategies..."
    if ! python scripts/comparison/01_compare_strategies.py; then
        print_error "Strategy comparison failed"
        return 1
    fi
    
    # Generate final report
    print_info "  Generating comprehensive report..."
    if ! Rscript scripts/comparison/02_generate_report.R; then
        print_warning "Report generation failed, but comparison data is available"
    fi
    
    mark_step_completed "$step_name"
    print_success "Strategy comparison and reporting completed"
}

# Function to display final summary
display_summary() {
    print_info "Pipeline execution completed!"
    print_info "=================================="
    
    echo ""
    echo "Results Summary:"
    echo "================"
    
    # Check which steps completed
    local completed_steps=()
    local failed_steps=()
    
    for step in "kmer_profiling" "variable_analysis" "coassembly_groups" "assemblies" "quality_assessment" "strategy_comparison"; do
        if check_step_completion "$step"; then
            completed_steps+=("$step")
        else
            failed_steps+=("$step")
        fi
    done
    
    echo "Completed steps (${#completed_steps[@]}): ${completed_steps[*]}"
    if [[ ${#failed_steps[@]} -gt 0 ]]; then
        echo "Failed/Skipped steps (${#failed_steps[@]}): ${failed_steps[*]}"
    fi
    
    echo ""
    echo "Key Output Files:"
    echo "=================="
    echo "• K-mer Analysis: results/similarity_analysis/"
    echo "• Variable Analysis: results/variable_analysis/"
    echo "• Assembly Results: results/assemblies/"
    echo "• Quality Assessment: results/quality_assessment/"
    echo "• Final Comparison: results/comparison/"
    echo "• Final Report: results/comparison/publication_figures/"
    
    # Try to display best strategy if available
    local recommendations_file="results/comparison/final_assembly_recommendations.txt"
    if [[ -f "$recommendations_file" ]]; then
        echo ""
        echo "Recommended Strategy:"
        echo "===================="
        grep "RECOMMENDED STRATEGY:" "$recommendations_file" | head -1
        grep "Overall Quality Score:" "$recommendations_file" | head -1
    fi
    
    echo ""
    echo "Next Steps:"
    echo "==========="
    echo "1. Review the final assembly recommendations"
    echo "2. Examine quality assessment results"
    echo "3. Use the recommended strategy for downstream analyses"
    echo "4. Consider running additional validation if needed"
}

# Function to clean up on exit
cleanup() {
    local exit_code=$?
    
    if [[ $exit_code -ne 0 ]]; then
        print_error "Pipeline failed with exit code $exit_code"
        print_info "Check log files in results/ for detailed error information"
    fi
    
    exit $exit_code
}

# Function to show usage
show_usage() {
    echo "Usage: $0 [config_file]"
    echo ""
    echo "Options:"
    echo "  config_file    Path to configuration file (default: config/config.yaml)"
    echo ""
    echo "Environment variables:"
    echo "  RESUME=true    Resume from last completed step"
    echo "  FORCE=true     Force re-run of all steps"
    echo ""
    echo "Examples:"
    echo "  $0                                    # Use default config"
    echo "  $0 my_config.yaml                    # Use custom config"
    echo "  RESUME=true $0                       # Resume interrupted run"
    echo "  FORCE=true $0                        # Force complete re-run"
}

# Main execution function
main() {
    # Set up error handling
    trap cleanup EXIT
    
    # Check for help request
    if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
        show_usage
        exit 0
    fi
    
    # Print welcome message
    print_info "Starting Viral Metagenomic Assembly Optimization Pipeline"
    print_info "Configuration: $CONFIG_FILE"
    print_info "Resume mode: $RESUME"
    print_info "Force mode: $FORCE"
    
    # Validate setup
    validate_config
    
    # Create results directory
    mkdir -p results
    
    # Execute pipeline steps
    run_kmer_profiling
    run_variable_analysis
    create_coassembly_groups
    run_assemblies
    run_quality_assessment
    run_strategy_comparison
    
    # Display final summary
    display_summary
    
    print_success "Viral metagenomic assembly optimization pipeline completed successfully!"
}

# Run main function with all arguments
main "$@"