# Anelloflow

`Anelloflow.sh` is a bash script designed for processing Oxford Nanopore Technologies (ONT) sequencing data, specifically targeting the assembly and identification of ICTV9 (https://ictv.global/report_9th/ssDNA/Anelloviridae) anellovirus species from metagenomic samples.

## Description

This script concatenates `fastq.gz` files from ONT sequencing runs, deconvolutes reads, and maps them against the ICTV9 anellovirus reference set. It filters out reads with a MAPQ score lower than 10, calculates mapping metrics, and reports the presence of anellovirus species and isolates in each sample at an abundance threshold of at least 5 reads.

### Command Line Usage

The script expects a directory of subdirectories as input. Each subdirectory should contain passing`fastq.gz` files for a specific sample/timepoint and be named after the sample/timepoint. 

```bash
./Anelloflow.sh /path/to/input_directory
```

## Dependencies

- `samtools`
- `minimap2`
- `seqkit`
- `ivar`

These dependencies must be installed and accessible in the system's `PATH`.

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourgithubusername/anelloflow.git
```

2. Ensure that the dependencies listed above are installed on your system.

3. Make the script executable:
```bash
chmod +x Anelloflow.sh
```

## Configuration

Before running the script, update the `USER DEFINED VARIABLES` section in the script to reflect your environment and project specifics:

- `HUMAN_REFERENCE:` Path to your human reference genome in `fasta.gz` format.
- `ANELLOVIRUS_REFERENCE:` Path to the Anellovirus reference set.
- `THREADS:` Number of CPU threads to use.

## Output

The script will generate the following files for each sample processed:

- A concatenated `fastq.gz` file of all passing reads.
- `BAM` files of reads mapped to the anellovirus reference set.
- A report of ICTV9 anellovirus species and isolates present in the sample.

## Notes

- This script is intended for use on Linux or macOS systems.
- Users should adjust the shebang line (`#!/opt/homebrew/bin/bash`) to match their system's environment if necessary.

## Citation

Anantharam, R., Duchen, D., Cox, A. L., Timp, W., Thomas, D. L., Clipman, S. J., & Kandathil, A. J. (2024). Long-Read Nanopore-Based Sequencing of Anelloviruses. Viruses, 16(5), 723. https://doi.org/10.3390/v16050723 

