# run_pipeline.py
import argparse
import os
import signal
import subprocess
import sys
from datetime import datetime
from pathlib import Path

child_proc = None  # Global variable to track the Snakemake process


def signal_handler(sig, frame):
    """Terminate Snakemake and its entire process group."""
    global child_proc
    print(
        f"\nReceived {signal.Signals(sig).name}, terminating Snakemake and subprocesses..."
    )
    if child_proc is not None:
        try:
            # Terminate the entire process group
            pgid = os.getpgid(child_proc.pid)
            os.killpg(pgid, sig)
            # Wait for process group to terminate
            child_proc.wait(timeout=15)
            print("Snakemake and subprocesses terminated.")
        except subprocess.TimeoutExpired:
            print("Forcing immediate termination with SIGKILL...")
            os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            print("Process group already terminated.")
        except Exception as e:
            print(f"Error terminating process group: {e}")
    sys.exit(1)


# Register signal handlers
signal.signal(signal.SIGTERM, signal_handler)
signal.signal(signal.SIGINT, signal_handler)


def run_snakemake(snakefile, profile, dry_run=False):
    """Run Snakemake with proper process group isolation."""
    global child_proc
    cmd = [
        "snakemake",
        "--snakefile",
        snakefile,
        "--profile",
        profile,
    ]
    if dry_run:
        cmd.append("-n")
        print(f"\n{' DRY RUN ':-^80}")
    else:
        print(f"\n{' RUNNING ':-^80}")

    print(f"Running command: {' '.join(cmd)}")

    try:
        # Start Snakemake in a new process group
        child_proc = subprocess.Popen(
            cmd,
            preexec_fn=os.setsid,  # Isolate process group
        )
        child_proc.wait()
        return child_proc.returncode == 0
    except Exception as e:
        print(f"\nException running {snakefile}: {e}")
        return False
    finally:
        child_proc = None  # Reset after completion


def main():
    parser = argparse.ArgumentParser(description="Run AlleleFlux pipeline.")
    parser.add_argument("--profile", required=True, help="Path to profile directory.")
    parser.add_argument("-n", "--dry-run", action="store_true", help="Dry run mode.")
    args = parser.parse_args()

    # Workflow 1: Preprocessing
    print("\n" + "=" * 80)
    print("Starting preprocessing workflow")
    print("=" * 80)
    if not run_snakemake("smk_workflow/snakefile.smk", args.profile, args.dry_run):
        print("\nPreprocessing failed! Aborting.")
        sys.exit(1)

    # Workflow 2: Allele analysis
    print("\n" + "=" * 80)
    print("Starting allele analysis workflow")
    print("=" * 80)
    if not run_snakemake("smk_workflow/snakefile2.smk", args.profile, args.dry_run):
        print("\nAllele analysis failed!")
        sys.exit(1)

    print("\nPipeline completed successfully!")


if __name__ == "__main__":
    main()
