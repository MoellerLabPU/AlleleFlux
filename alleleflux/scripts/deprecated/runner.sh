#!/bin/bash

# AlleleFlux Sequential Workflow Runner
# =====================================
# 
# This script sequentially executes step1.smk and step2.smk with proper error handling,
# signal management, and SLURM job cleanup.

# Bash strict mode: exit on error, undefined variables, and pipe failures
set -euo pipefail

# Filename of the script (for usage messages)
SCRIPT_NAME=$(basename "$0")

# --- Path Configuration: Make the script location-aware ---
# Get the absolute directory where this script is located. This is the crucial
# step that makes the script portable and runnable from any location.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

# Define absolute paths to the bundled workflow files.
WORKFLOW_DIR="$SCRIPT_DIR/smk_workflow"
SNAKEFILE_STEP1="$WORKFLOW_DIR/step1.smk"
SNAKEFILE_STEP2="$WORKFLOW_DIR/step2.smk"
DEFAULT_CONFIG_FILE="$WORKFLOW_DIR/config.yml"

# Default configuration values - can be overridden by command line arguments
CONFIG_FILE="$DEFAULT_CONFIG_FILE"      # Path to Snakemake config file
CORES=1                                 # Number of CPU cores to use
PROFILE=""                              # Snakemake profile directory (empty = local execution)
DRY_RUN=false                          # Whether to perform a dry run (no actual execution)
STEP1_ONLY=false                       # Run only Step 1 if true
STEP2_ONLY=false                       # Run only Step 2 if true

# ANSI color codes for colored terminal output
RED='\033[0;31m'      
GREEN='\033[0;32m'   
YELLOW='\033[1;33m'   
BLUE='\033[0;34m'     
NC='\033[0m'          # No Color - reset to default

# Logging functions with color-coded output
# All log messages are sent to stderr (&2) to avoid interfering with script output

# Blue text for info messages
log_info() {
    echo -e "${BLUE}[INFO]${NC} $*" >&2
}

 # Green text for success messages
log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $*" >&2
}

# Yellow text for warnings
log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $*" >&2
}

# Red text for errors
log_error() {
    echo -e "${RED}[ERROR]${NC} $*" >&2
}

# Function to display usage
usage() {
    cat << EOF
Usage: $SCRIPT_NAME [OPTIONS]

This script runs AlleleFlux workflows sequentially:
1. Step 1: Sample profiling and eligibility table generation
2. Step 2: Allele analysis and scoring (waits for Step 1 completion)

Options:
    -c, --config FILE       Configuration file (default: $CONFIG_FILE)
    -j, --cores N           Number of cores to use (default: $CORES)
    -p, --profile DIR       Snakemake profile directory for cluster execution
    -n, --dry-run           Perform a dry run without executing jobs
    -1, --step1-only        Run only Step 1
    -2, --step2-only        Run only Step 2 (assumes Step 1 completed)
    -h, --help              Show this help message

Examples:
    # Run complete workflow with 8 cores
    $SCRIPT_NAME --config smk_workflow/config.yml --cores 8

    # Run with SLURM profile on cluster
    $SCRIPT_NAME --profile ~/smk_workflow/profile --cores 100

    # Dry run to check workflow
    $SCRIPT_NAME -n

    # Run only Step 1
    $SCRIPT_NAME --step1-only --cores 16

    # Run only Step 2 (Step 1 must be completed)
    $SCRIPT_NAME --step2-only --cores 8

EOF
}


# Parse command line arguments using a while loop
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--config)
            CONFIG_FILE="$2" 
            shift 2           
            ;;
        -j|--cores)
            CORES="$2"        
            shift 2
            ;;
        -p|--profile)
            PROFILE="$2"      
            shift 2
            ;;
        -n|--dry-run)
            DRY_RUN=true   
            shift
            ;;
        -1|--step1-only)
            STEP1_ONLY=true
            shift
            ;;
        -2|--step2-only)
            STEP2_ONLY=true
            shift
            ;;
        -h|--help)
            usage            # Show help and exit
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Validate that conflicting options aren't specified together
if [[ "$STEP1_ONLY" == true && "$STEP2_ONLY" == true ]]; then
    log_error "Cannot specify both --step1-only and --step2-only"
    exit 1
fi

# Validate that cores is a positive integer
if ! [[ "$CORES" =~ ^[0-9]+$ ]] || [[ "$CORES" -lt 1 ]]; then
    log_error "Cores must be a positive integer, got: $CORES"
    exit 1
fi

# Function to validate that all required prerequisites are available
validate_prerequisites() {
    # Check if snakemake command is installed and available in PATH
    if ! command -v snakemake &> /dev/null; then
        log_error "snakemake command not found. Please install Snakemake."
        exit 1
    fi
    
    # Check if the specified config file exists
    if [[ ! -f "$CONFIG_FILE" ]]; then
        log_error "Configuration file not found: $CONFIG_FILE"
        exit 1
    fi
    
    # Check if we're running from the correct directory (AlleleFlux root)
    if [[ ! -d "$WORKFLOW_DIR" ]]; then
        log_error "Packaged smk_workflow directory not found at: $WORKFLOW_DIR"
        log_error "This may indicate a broken installation."
        exit 1
    fi
    
    log_info "Prerequisites validated successfully"
}

# Function to clean up resources when the script is interrupted or exits
# This ensures no orphaned processes or SLURM jobs are left running
cleanup() {
    local exit_code=$?  # Capture the exit code from the interrupted process
    log_warning "Caught interrupt signal. Cleaning up..."
    
     # Kill any running snakemake processes that we started and their child processes
    if [[ -n "${SNAKEMAKE_PID:-}" ]]; then
        log_info "Terminating Snakemake process $SNAKEMAKE_PID and its children..."
        
        # Kill the entire process group to catch child processes (like Python multiprocessing workers)
        # First try graceful termination of the process group
        kill -TERM -"$SNAKEMAKE_PID" 2>/dev/null || true
        sleep 5  # Give more time for Python multiprocessing cleanup
        
        # Force kill the entire process group if still running
        kill -KILL -"$SNAKEMAKE_PID" 2>/dev/null || true
        
        # Also try to kill any remaining alleleflux processes that might be orphaned
        log_info "Cleaning up any remaining AlleleFlux processes..."
        pkill -f "alleleflux-" 2>/dev/null || true

    fi
    
    # If SLURM commands are available, cancel all jobs for current user
    # This prevents orphaned jobs from continuing to run on the cluster
    if command -v scancel &> /dev/null && command -v squeue &> /dev/null; then
        log_info "Cancelling SLURM jobs for user $USER..."
        local job_ids
        # Get all job IDs for the current user
        job_ids=$(squeue -u "$USER" -h -o "%i" 2>/dev/null || true)
        if [[ -n "$job_ids" ]]; then
            # Cancel each job ID
            echo "$job_ids" | xargs -r scancel 2>/dev/null || true
            log_info "Cancelled $(echo "$job_ids" | wc -w) SLURM jobs"
        else
            log_info "No SLURM jobs found for user $USER"
        fi
    fi
    
    log_info "Cleanup completed."
    exit $exit_code  # Exit with the original exit code
}

# Set up signal traps to catch interruptions (Ctrl+C) and script termination
# SIGINT = Ctrl+C, SIGTERM = termination signal, EXIT = any script exit
trap cleanup SIGINT SIGTERM EXIT

# Function to safely unlock Snakemake directory if needed
# This handles stale lock files that can occur after kill signals or power loss
safe_snakemake_unlock() {
    local workflow_dir="${1:-$(pwd)}"  # Use provided directory or current working directory
    local snakefile="$2"
    local lock_dir="$workflow_dir/.snakemake/locks"
    
    log_info "Checking for Snakemake locks in directory: $workflow_dir"
    
    # Check if lock directory exists (if not, no unlock needed)
    if [[ ! -d "$lock_dir" ]]; then
        log_info "No Snakemake lock directory found. No unlock needed."
        return 0  # No locks exist, nothing to do
    fi
    
    # Check if any lock files exist
    if ! find "$lock_dir" -name "*.lock" -type f 2>/dev/null | grep -q .; then
        log_info "No lock files found in lock directory. No unlock needed."
        return 0  # No lock files found, nothing to do
    fi
    
    log_info "Detected Snakemake lock files. Checking for active Snakemake processes..."
    
    # Check for active Snakemake processes using multiple methods for robustness
    local snakemake_running="false"
    
    # Method 1: Check if SNAKEMAKE_PID is set and that process is still running
    # This won't work as we are running this function before starting any Snakemake jobs
    # if [[ -n "${SNAKEMAKE_PID:-}" ]]; then
    #     log_info "Checking if tracked Snakemake process $SNAKEMAKE_PID is still running..."
    #     if ps -p "$SNAKEMAKE_PID" >/dev/null 2>&1; then
    #         log_warning "Snakemake process $SNAKEMAKE_PID is still running. Will not unlock."
    #         snakemake_running="true"
    #     else
    #         log_info "SNAKEMAKE_PID $SNAKEMAKE_PID is no longer running (stale PID)"
    #         unset SNAKEMAKE_PID  # Clear stale PID
    #     fi
    # fi
    
    # Method 2: Check for any Snakemake processes by name with more precise pattern
    if [[ "$snakemake_running" == "false" ]]; then
        log_info "Searching for active Snakemake processes using pattern matching..."
        local snakemake_pids
        # Use more precise pattern to avoid catching unrelated Snakemake jobs
        # This is an alternative to ps -u suppal -o pid,cmd |   grep -E "snakemake.*--snakefile.*smk_workflow" | grep -v grep
        snakemake_pids=$(pgrep -u "$USER" -f "snakemake.*--snakefile.*smk_workflow" 2>/dev/null || true)
        if [[ -n "$snakemake_pids" ]]; then
            log_warning "Found active Snakemake processes with smk_workflow pattern (PIDs: $snakemake_pids). Will not unlock."
            snakemake_running="true"
        else
            log_info "No active smk_workflow Snakemake processes found."
        fi
    fi
    
    # Method 3: Check for processes with open files in the lock directory
    if [[ "$snakemake_running" == "false" ]] && command -v lsof >/dev/null 2>&1; then
        log_info "Checking for processes with open files in lock directory..."
        local lock_processes
        lock_processes=$(lsof +D "$lock_dir" 2>/dev/null | awk 'NR>1 {print $2}' | sort -u || true)
        if [[ -n "$lock_processes" ]]; then
            log_warning "Found processes with open files in lock directory (PIDs: $lock_processes). Will not unlock."
            snakemake_running="true"
        else
            log_info "No processes found with open files in lock directory."
        fi
    fi
    
    # If no active Snakemake processes detected, proceed with safe unlock
    if [[ "$snakemake_running" == "false" ]]; then
        log_info "No active Snakemake processes detected. Attempting safe unlock..."
        
        # Method 1: Try snakemake --unlock first (safest approach)
        if command -v snakemake >/dev/null 2>&1; then
            log_info "Running: snakemake --unlock --snakefile "$snakefile" --cores 1 --directory $workflow_dir"
            if snakemake --unlock --snakefile "$snakefile" --cores 1 --directory "$workflow_dir" 2>/dev/null; then
                log_success "Successfully unlocked Snakemake directory using --unlock"
                return 0
            else
                log_warning "snakemake --unlock failed. Attempting manual lock cleanup..."
            fi
        else
            log_warning "snakemake command not available. Attempting manual lock cleanup..."
        fi
        
        # Method 2: Manual cleanup as fallback (only if snakemake --unlock fails)
        # log_info "Manually removing lock files from $lock_dir"
        # local lock_count
        # lock_count=$(find "$lock_dir" -name "*.lock" -type f 2>/dev/null | wc -l)
        
        # if [[ "$lock_count" -gt 0 ]]; then
        #     # List the lock files before removal for transparency
        #     log_info "Found $lock_count lock file(s):"
        #     find "$lock_dir" -name "*.lock" -type f 2>/dev/null | while read -r lock_file; do
        #         log_info "  - $(basename "$lock_file")"
        #     done
            
        #     # Remove the lock files
        #     if find "$lock_dir" -name "*.lock" -type f -delete 2>/dev/null; then
        #         log_success "Successfully removed $lock_count stale lock file(s)"
        #     else
        #         log_error "Failed to remove some lock files. Manual intervention may be required."
        #         return 1
        #     fi
        # fi
        
        log_success "Snakemake unlock completed successfully"
        return 0
        
    else
        log_warning "Skipping unlock due to active Snakemake processes"
        log_warning "If you're certain no Snakemake is running, manually run: snakemake --unlock"
        return 1
    fi
}

# Function to run a single Snakemake workflow step
# Parameters: step_name, snakefile_path
run_snakemake_steps() {
    local step_name="$1"    # Human-readable name for logging
    local snakefile="$2"    # Path to the Snakemake file
    
    # Log the start of this step with key parameters
    log_info "Starting $step_name..."
    log_info "Snakefile: $snakefile"
    log_info "Config: $CONFIG_FILE"
    log_info "Cores: $CORES"
    
    if [[ -n "$PROFILE" ]]; then
        log_info "Profile: $PROFILE"
    fi
    
    # Validate that the Snakefile exists before trying to run it
    if [[ ! -f "$snakefile" ]]; then
        log_error "Snakefile not found: $snakefile"
        exit 1
    fi
    
    # Safely unlock Snakemake directory if there are stale locks
    if ! safe_snakemake_unlock "$(pwd)" "$snakefile"; then
        log_error "Failed to safely unlock Snakemake directory. Aborting $step_name."
        exit 1
    fi
    
    # Build the snakemake command as an array for proper argument handling
    local cmd=(snakemake)
    cmd+=("--snakefile" "$snakefile")           # Specify which Snakefile to use
    cmd+=("--configfile" "$CONFIG_FILE")        # Specify config file
    # cmd+=("--cores" "$CORES")                   # Number of cores to use
    cmd+=("--rerun-incomplete")                 # Rerun any incomplete jobs from previous runs
    
    # Add profile if specified (for cluster execution)
    if [[ -n "$PROFILE" ]]; then
        cmd+=("--profile" "$PROFILE")
    fi
    
    # Conditionally add --cores.
    # Add it if No profile is used (local execution needs it).
    if [[ -z "$PROFILE" ]]; then
        cmd+=("--cores" "$CORES")
    fi
    
    # Add dry run flag if requested (test mode - no actual execution)
    if [[ "$DRY_RUN" == true ]]; then
        cmd+=("--dry-run")
        log_warning "DRY RUN MODE - No jobs will be executed"
    fi
    
    # Display the full command being executed
    log_info "Command: ${cmd[*]}"
    echo "========================================"
    
    # Execute the command and track its execution
    local start_time=$(date +%s)  # Record start time for duration calculation
    SNAKEMAKE_PID=""               # Initialize PID variable
    
    # Run the command in background and capture its PID
    # Use setsid to create a new process group for proper cleanup
    if setsid "${cmd[@]}" & then
        SNAKEMAKE_PID=$!           # Store the process ID for potential cleanup
        export SNAKEMAKE_PID       # Export PID so it's available to other functions
        log_info "Started Snakemake process with PID: $SNAKEMAKE_PID"
        wait $SNAKEMAKE_PID        # Wait for the process to complete
        local exit_code=$?         # Capture the exit code
        SNAKEMAKE_PID=""           # Clear PID since process is done
        unset SNAKEMAKE_PID        # Also unset the exported variable
        
        # Calculate execution duration
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        
        # Check if the step completed successfully
        if [[ $exit_code -eq 0 ]]; then
            log_success "$step_name completed successfully in ${duration}s"
            echo "========================================"
            return 0
        else
            log_error "$step_name failed with exit code $exit_code after ${duration}s"
            return $exit_code
        fi
    else
        log_error "Failed to start $step_name"
        return 1
    fi
}

# Main execution function - orchestrates the entire workflow
main() {
    local overall_start_time=$(date +%s)  # Record start time for total duration
    
    # Log startup information
    log_info "AlleleFlux Sequential Workflow Runner Starting"
    log_info "Version: ${ALLELEFLUX_VERSION:-unknown}" # Use version from env var
    log_info "Packaged content is at: $SCRIPT_DIR"
    log_info "Timestamp: $(date)"
    log_info "Working directory: $(pwd)"
    log_info "User: $USER"
    
    # Validate that all prerequisites are met before starting
    validate_prerequisites
    
    # Execute workflows based on the options specified
    if [[ "$STEP2_ONLY" == true ]]; then
        # User requested Step 2 only - assumes Step 1 outputs already exist
        log_info "Running Step 2 only (assuming Step 1 is completed)"
        run_snakemake_steps "Step 2: Allele Analysis and Scoring" "$SNAKEFILE_STEP2"
        
    elif [[ "$STEP1_ONLY" == true ]]; then
        # User requested Step 1 only - useful for testing or partial runs
        log_info "Running Step 1 only"
        run_snakemake_steps "Step 1: Sample Profiling and Eligibility" "$SNAKEFILE_STEP1"
        
    else
        # Default: run complete workflow sequentially
        log_info "Running complete workflow (Step 1 → Step 2)"
        
        # Run Step 1: Sample profiling and eligibility table generation
        log_info "========== STEP 1 =========="
        run_snakemake_steps "Step 1: Sample Profiling and Eligibility" "$SNAKEFILE_STEP1"
        
        # Only proceed to Step 2 if Step 1 succeeded (run_snakemake_steps would exit on failure)
        log_info "Step 1 completed successfully. Proceeding to Step 2..."
        sleep 2  # Brief pause for clarity in logs
        
        # Run Step 2: Allele analysis and scoring (depends on Step 1 outputs)
        log_info "========== STEP 2 =========="
        run_snakemake_steps "Step 2: Allele Analysis and Scoring" "$SNAKEFILE_STEP2"
    fi
    
    # Calculate and display total execution time
    local overall_end_time=$(date +%s)
    local total_duration=$((overall_end_time - overall_start_time))
    
    # Log successful completion
    log_success "AlleleFlux workflow completed successfully!"
    log_info "Total execution time: ${total_duration}s"
    log_info "Finished at: $(date)"
}

# Execute the main function with all command line arguments
# "$@" passes all arguments that were given to the script
main "$@"
