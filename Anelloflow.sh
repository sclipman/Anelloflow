#!/opt/homebrew/bin/bash

# ------------------------------------------------------------------------------------------------------------------------------------
# Anelloflow.sh : Anellovirus Virome Assembly Pipline
#
# Author(s): Steven J. Clipman (sclipman@jhmi.edu)
#
# Version: 1.0
#
# Description: This script processes passing fastq.gz files from ONT sequencing. 
# Command line input is a directory of subdirectories containing passing fastq.gz for each sample/timepoint.
# The subdirectories should be named with the sample name/timepoint. The script will concatenates fastq.gz files, deconvolute reads,
# and map them to the ICTV9 anellovirus reference. The script then filters reads with MAPQ < 10, calculates mapping metrics, and
# reports the presence of anellovirus species and isolates in each sample at an abundance of at least 5 reads.
#
# Dependencies: samtools, minimap2, seqkit, ivar
#
# Notes: This script is designed to be run on Linux or Mac OS. The user should update the shebang line to the path of their bash.
# The user should also updated variables in the "USER DEFINED VARIABLES" section below.
# ------------------------------------------------------------------------------------------------------------------------------------


# USER DEFINED VARIABLES
# -------------------------------------------------------------------------------------

# Human reference
HUMAN_REFERENCE=""

# Anellovirus reference
ANELLOVIRUS_REFERENCE=""

# Number of CPU Threads
THREADS=20


# Functions
# -------------------------------------------------------------------------------------

# Define function to check dependencies
check_dependencies() {
    local missing_counter=0
    for dependency in samtools minimap2 seqkit ivar; do
        if ! command -v $dependency &> /dev/null; then
            echo_red "Error: $dependency is not installed."
            ((missing_counter++))
        fi
    done

    if [ $missing_counter -ne 0 ]; then
        echo_red "Error: $missing_counter dependencies are missing. Install them before running this script."
        exit 1
    else
        echo_green "All dependencies are installed."
    fi

}

# Define function to process species presence
process_species_presence() {
    local INPUT_BAM="$1"
    local OUTPUT_FILE="$2"
    local REFS_FILE="$3"
        echo ""
        echo "Anellovirus species and isolates detected in $INPUT_BAM:"

    # Initialize an associative array to hold species and isolate names
    declare -A species_map
    species_map["AB008394"]="Torque teno virus 1 - TA278"
    species_map["AB017613"]="Torque teno virus 16 - TUS01"
    species_map["AB025946"]="Torque teno virus 19 - SANBAN"
    species_map["AB026929"]="Torque teno mini virus 6 - CBD203"
    species_map["AB026931"]="Torque teno mini virus 1 - CBD279"
    species_map["AB028668"]="Torque teno virus 15 - TJN01"
    species_map["AB037926"]="Torque teno virus 14 - CH65-1"
    species_map["AB038621"]="Torque teno virus 29 - yonKC009"
    species_map["AB038627"]="Torque teno mini virus 7 - CLC156"
    species_map["AB038629"]="Torque teno mini virus 2 - NLC023"
    species_map["AB038630"]="Torque teno mini virus 3 - NLC026"
    species_map["AB038631"]="Torque teno mini virus 9 - NLC030"
    species_map["AB041957"]="Torque teno virus 4 - Pt-TTV6"
    species_map["AB041958"]="Torque teno virus 26 - Mf-TTV3"
    species_map["AB041959"]="Torque teno virus 25 - Mf-TTV9"
    species_map["AB041960"]="Torque teno tamarin virus - So-TTV2"
    species_map["AB041961"]="Torque teno douroucouli virus - At-TTV3"
    species_map["AB041962"]="Torque teno mini virus 5 - TGP96"
    species_map["AB041963"]="Torque teno mini virus 4 - Pt-TTV8-II"
    species_map["AB049607"]="Torque teno virus 23 - CH65-2"
    species_map["AB049608"]="Torque teno virus 2 - CH71"
    species_map["AB054647"]="Torque teno virus 8 - Kt-08F"
    species_map["AB057358"]="Torque teno tupaia virus - Tbc-TTV14"
    species_map["AB060594"]="Torque teno virus 20 - SAa-10"
    species_map["AB060597"]="Torque teno virus 24 - SAa-01"
    species_map["AB064595"]="Torque teno virus 27 - CT23F"
    species_map["AB064598"]="Torque teno virus 28 - CT43F"
    species_map["AB064605"]="Torque teno virus 12 - CT44F"
    species_map["AB064607"]="Torque teno virus 10 - JT34F"
    species_map["AB076001"]="Torque teno sus virus 1 - Sd-TTV31"
    species_map["AB076002"]="Torque teno canis virus - Cf-TTV10"
    species_map["AB076003"]="Torque teno felis virus - Fc-TTV4"
    species_map["AB290918"]="Torque teno midi virus 1 - MD1-073"
    species_map["AB290919"]="Torque teno midi virus 2 - MD2-013"
    species_map["AB303552"]="Torque teno midi virus - MDJHem2"
    species_map["AB303553"]="Torque teno midi virus - MDJHem3-1"
    species_map["AB303554"]="Torque teno midi virus - MDJHem3-2"
    species_map["AB303555"]="Torque teno midi virus - MDJHem5"
    species_map["AB303559"]="Torque teno midi virus - MDJN2"
    species_map["AB303560"]="Torque teno midi virus - MDJN14"
    species_map["AB303561"]="Torque teno midi virus - MDJN47"
    species_map["AB303562"]="Torque teno midi virus - MDJN51"
    species_map["AB303564"]="Torque teno midi virus - MDJN69"
    species_map["AB303566"]="Torque teno midi virus - MDJN97"
    species_map["AB449062"]="Torque teno midi virus - Pt-TTMDV210"
    species_map["AF261761"]="Torque teno virus 7 - PMV"
    species_map["AF291073"]="Torque teno mini virus 8 - PB4TL"
    species_map["AF345523"]="Torque teno virus 5 - TCHN-C1"
    species_map["AF345524"]="Torque teno virus 11 - TCHN-D1"
    species_map["AF345526"]="Torque teno virus 13 - TCHN-A"
    species_map["AF348409"]="Torque teno virus 21 - TCHN-B"
    species_map["AF435014"]="Torque teno virus 6 - KAV"
    species_map["AX025718"]="Torque teno virus 18 - SENV-C"
    species_map["AX025830"]="Torque teno virus 17 - SENV-G"
    species_map["AX174942"]="Torque teno virus 22 - svi-1"
    species_map["AY666122"]="Torque teno virus 3 - HEL32"
    species_map["AY823990"]="Torque teno sus virus 2 - 1p"
    species_map["AY823991"]="Torque teno sus virus - 2p"
    species_map["DQ187006"]="Torque teno virus 9 - BM1C-18"
    species_map["EF538875"]="Torque teno midi virus - 2PoSMA"
    species_map["EF538876"]="Torque teno midi virus - 6PoSMA"
    species_map["EF538877"]="Torque teno felis virus - PRA1"
    species_map["EF538880"]="Torque teno mini virus - LIL-y1"
    species_map["EF538881"]="Torque teno mini virus - LIL-y2"
    species_map["EF538882"]="Torque teno mini virus - LIL-y3"
    species_map["FJ459582"]="Torque teno zalophus virus – ZcAV"

# Initialize an associative array to hold references with at least 5 reads
    declare -A refs_with_min_reads
    while IFS= read -r ref; do
        # Extract just the accession number before the underscore
        ref=$(echo "$ref" | cut -d '_' -f 1)
        refs_with_min_reads["$ref"]=1
    done < "$REFS_FILE"

    # Append the detected species and isolates to both terminal and output file
    samtools idxstats "$INPUT_BAM" | while IFS=$'\t' read -r ref count _ _; do
        # Extract just the accession number before the underscore for comparison
        ref=$(echo "$ref" | cut -d '_' -f 1)
        if [[ ${refs_with_min_reads["$ref"]} ]]; then
            if [[ -n "${species_map["$ref"]}" ]]; then
                echo "$ref    ${species_map["$ref"]}"
            else
                echo "Unknown or unmapped for accession $ref"
            fi
        fi
    done | sort | uniq | tee -a "$OUTPUT_FILE"

}

# Define function to generate consensus sequences with ivar
generate_consensus() {
    local INPUT_BAM="$1" # BAM file filtered by MAPQ
    local REFS_FILE="$2" # File containing reference IDs for detected viruses
    local SAMPLE_DIR="$3" # Directory for storing results
    local SAMPLE_NAME="$4" # Name of the sample being processed
    local REFERENCE="$5" # Anellovirus reference genome
    
    mkdir -p "${SAMPLE_DIR}/consensus_sequences"

    # Ensure the BAM file is indexed for samtools view to work properly
    samtools index "$INPUT_BAM"

    # Iterate over each reference (virus) identified in REFS_FILE
    while IFS= read -r ref; do
        local REF_ACCESSION="${ref%%_*}" # Extract accession number before the underscore if present
        local CONSENSUS_OUTPUT="${SAMPLE_DIR}/consensus_sequences/${SAMPLE_NAME}_${REF_ACCESSION}_consensus.fa"
        
        echo_blue "Generating consensus for $REF_ACCESSION..."

        # Filter BAM for reads specific to this reference sequence
        samtools view -b "$INPUT_BAM" "$REF_ACCESSION" > "${SAMPLE_DIR}/${SAMPLE_NAME}_${REF_ACCESSION}_filtered.bam"
        
        # Check if the filtered BAM file is non-empty before proceeding
        if [ -s "${SAMPLE_DIR}/${SAMPLE_NAME}_${REF_ACCESSION}_filtered.bam" ]; then
            samtools mpileup -A -d 0 -Q 0 --reference "$REFERENCE" "${SAMPLE_DIR}/${SAMPLE_NAME}_${REF_ACCESSION}_filtered.bam" | \
            ivar consensus -p "${CONSENSUS_OUTPUT%.fa}" -n N -q 10 -m 1

            if [ -f "${CONSENSUS_OUTPUT}" ]; then
                echo_green "Consensus sequence generated for $REF_ACCESSION: $CONSENSUS_OUTPUT"
            else
                echo_red "Failed to generate consensus for $REF_ACCESSION."
            fi
        else
            echo_yellow "No reads for $REF_ACCESSION meeting criteria; skipping consensus generation."
        fi

        # Cleanup to save space
        rm "${SAMPLE_DIR}/${SAMPLE_NAME}_${REF_ACCESSION}_filtered.bam" "${SAMPLE_DIR}/${SAMPLE_NAME}_${REF_ACCESSION}_filtered.bam.bai"

    done < "$REFS_FILE"
}

# Define terminal print colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print in red
echo_red() {
    echo -e "${RED}$1${NC}"
}

# Function to print in green
echo_green() {
    echo -e "${GREEN}$1${NC}"
}

# Function to print in yellow
echo_yellow() {
    echo -e "${YELLOW}$1${NC}"
}

# Function to print in blue
echo_blue() {
    echo -e "${BLUE}$1${NC}"
}

# Check for the correct usage and presence of an argument
if [ "$#" -ne 1 ]; then
    echo_red "Error: Incorrect usage."
    echo "Usage: $0 input_directory"
    exit 1
fi


# Main Pipeline
# -------------------------------------------------------------------------------------

# Check for dependencies
echo_blue "Checking for dependencies..."
check_dependencies

# Assign the first argument to INPUT_DIRECTORY
INPUT_DIRECTORY="$1"

# Reference indices
HUMAN_REFERENCE_INDEX="${HUMAN_REFERENCE}.mmi"
ANELLOVIRUS_REFERENCE_INDEX="${ANELLOVIRUS_REFERENCE}.mmi"
ANELLOVIRUS_REFERENCE_FAI="${ANELLOVIRUS_REFERENCE}.fai"

# Index references if needed
echo ""
echo_blue "Checking for reference indices"
index_and_print() {
    local REFERENCE_INDEX=$1
    local REFERENCE=$2
    if [ ! -f "$REFERENCE_INDEX" ]; then
        echo_yellow "Minimap index not found. Indexing $REFERENCE"
        minimap2 -d "$REFERENCE_INDEX" "$REFERENCE"
        echo_green "Indexing of $REFERENCE completed."
    else
        echo_green "Minimap index for $REFERENCE is available."
    fi
}

index_and_print "$HUMAN_REFERENCE_INDEX" "$HUMAN_REFERENCE" 
index_and_print "$ANELLOVIRUS_REFERENCE_INDEX" "$ANELLOVIRUS_REFERENCE"


# Process each sample subdirectory
echo ""
echo_blue "Beginning sample processing..."
for SAMPLE_DIR in "${INPUT_DIRECTORY}"/*/; do

    # Trim trailing slash for consistency
    SAMPLE_DIR=${SAMPLE_DIR%/}
    if [ -d "$SAMPLE_DIR" ] && [ "$(ls -A "$SAMPLE_DIR")" ]; then
        echo ""
        echo_blue "Processing directory: $SAMPLE_DIR"
        FULL_SAMPLE_NAME=$(basename "$SAMPLE_DIR")
        
        FASTQ_PATTERN="${SAMPLE_DIR}"/*.fastq.gz
        shopt -s nullglob
        FASTQ_FILES=($FASTQ_PATTERN)
        shopt -u nullglob
        
        COMBINED_FASTQ="${SAMPLE_DIR}/${FULL_SAMPLE_NAME}_all_passing_reads.fastq.gz"
        if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
            echo_yellow "No .fastq.gz files found in $SAMPLE_DIR"
            continue
        else
            echo "Found ${#FASTQ_FILES[@]} .fastq.gz files in $SAMPLE_DIR"
            if [ -f "$COMBINED_FASTQ" ]; then
                echo_green "Concatenated fastq file already exists for $FULL_SAMPLE_NAME, skipping concatenation."
            else
                echo ""
                echo_blue "Concatenating files into $COMBINED_FASTQ"
                gunzip -c "${FASTQ_FILES[@]}" | gzip > "$COMBINED_FASTQ"

                #Remove individual fastq files
                echo_blue "Removing individual fastq files"
                rm "${FASTQ_FILES[@]}"
            fi
        fi

        # Separate reads into short and long reads after concatenation
        echo_blue "Separating reads into short (< 1000 bp) and long (> 1000 bp) reads for $FULL_SAMPLE_NAME"
        seqkit seq -m 1000 "$COMBINED_FASTQ" | gzip > "${SAMPLE_DIR}/${FULL_SAMPLE_NAME}_all_passing_reads.long.fastq.gz"
        seqkit seq -M 999 "$COMBINED_FASTQ" | gzip > "${SAMPLE_DIR}/${FULL_SAMPLE_NAME}_all_passing_reads.short.fastq.gz"

        # Perform host deconvolution for long and short reads
        for READ_TYPE in short long; do
            DECON_OUTPUT="${SAMPLE_DIR}/${FULL_SAMPLE_NAME}_hostdecon.${READ_TYPE}.fastq.gz"
            
            # Check if host deconvoluted reads already exist
            if [ -f "$DECON_OUTPUT" ]; then
                echo_green "Host deconvoluted ${READ_TYPE} reads already exist for $FULL_SAMPLE_NAME, skipping."
                continue # Skip to the next loop iteration if file exists
            fi
            
            echo ""
            echo_blue "Host deconvoluting ${READ_TYPE} reads for $FULL_SAMPLE_NAME."
            
            # Use the appropriate minimap2 preset for each read type
            if [[ "$READ_TYPE" == "short" ]]; then
                MINIMAP2_PRESET="sr"
            else
                MINIMAP2_PRESET="map-ont"
            fi
            
            minimap2 -ax ${MINIMAP2_PRESET} -t "$THREADS" --secondary=no "$HUMAN_REFERENCE_INDEX" "${SAMPLE_DIR}/${FULL_SAMPLE_NAME}_all_passing_reads.${READ_TYPE}.fastq.gz" \
            | samtools view -b -f 4 - \
            | samtools sort -n - \
            | samtools fastq - \
            | gzip > "$DECON_OUTPUT"
        done


        # Map host deconvoluted short and long reads to the anellovirus reference set
        echo ""
        echo_blue "Mapping deconvoluted reads for $SAMPLE_DIR to anellovirus references..."
        
        # Initialize an array to store BAM file names for later combination
        declare -a BAM_FILES=()
        
        for READ_TYPE in short long; do
            ANELLOMAP_OUTPUT="${SAMPLE_DIR}/${FULL_SAMPLE_NAME}_anellovirus_${READ_TYPE}.bam"
            if [[ "$READ_TYPE" == "short" ]]; then
                MINIMAP2_PRESET="sr"
            else
                MINIMAP2_PRESET="map-ont"
            fi
            
            minimap2 -ax ${MINIMAP2_PRESET} -t $THREADS --secondary=no "$ANELLOVIRUS_REFERENCE_INDEX" "${SAMPLE_DIR}/${FULL_SAMPLE_NAME}_hostdecon.${READ_TYPE}.fastq.gz" \
            | samtools view -b -o "$ANELLOMAP_OUTPUT"
            
            BAM_FILES+=("$ANELLOMAP_OUTPUT") # Add output BAM to array for later combination
            
            echo_green "Mapping complete for $FULL_SAMPLE_NAME ($READ_TYPE reads)"
        done

        # Combine short and long read BAM files
        COMBINED_BAM="${SAMPLE_DIR}/${FULL_SAMPLE_NAME}_anellovirus_combined.bam"
        samtools merge "$COMBINED_BAM" "${BAM_FILES[@]}"
        echo_green "Combined BAM file created for $FULL_SAMPLE_NAME."

        # Filter out reads with MAPQ < 10 from the combined BAM file
        ANELLOMAP_FILTERED="${SAMPLE_DIR}/${FULL_SAMPLE_NAME}_anellovirus_filtered.bam"
        echo ""
        echo_blue "Filtering reads with MAPQ < 10 for $FULL_SAMPLE_NAME"
        samtools view -b -q 10 "$COMBINED_BAM" > "$ANELLOMAP_FILTERED"
        echo_green "Filtering complete for $FULL_SAMPLE_NAME."

        # Sorting and Indexing the filtered BAM file
        SORTED_ANELLOMAP_FILTERED="${SAMPLE_DIR}/${FULL_SAMPLE_NAME}_anellovirus_filtered_sorted.bam"
        samtools sort "$ANELLOMAP_FILTERED" -o "$SORTED_ANELLOMAP_FILTERED"
        samtools index "$SORTED_ANELLOMAP_FILTERED"
        echo_green "Filtering complete for $FULL_SAMPLE_NAME"

        # Calculate depth to later filter based on 5X coverage
            # Count reads per reference
            samtools idxstats "$SORTED_ANELLOMAP_FILTERED" > "${SAMPLE_DIR}/read_counts_per_ref.txt"

            # Extract references with >= 5 reads
            awk '$3 >= 5 {print $1}' "${SAMPLE_DIR}/read_counts_per_ref.txt" > "${SAMPLE_DIR}/refs_with_min_5_reads.txt"

    else
        echo_red "Error: Directory $SAMPLE_DIR does not exist or is empty."
    fi

    # Calculate mapping metrics
    echo ""
    echo_blue "Getting anellovirus mapping metrics for $FULL_SAMPLE_NAME"
    
    # Define the input BAM file name
    INPUT_BAM="$SORTED_ANELLOMAP_FILTERED"

    # Extract mapping qualities using samtools and store them in an array
    MAP_QUALS=($(samtools view "$INPUT_BAM" | awk '{print $5}' | sort -n))

    if [ ${#MAP_QUALS[@]} -eq 0 ]; then
        echo "No mapped reads."
        exit 1
    fi

    # Mean calculation
    MEAN=$(printf "%s\n" "${MAP_QUALS[@]}" | awk '{sum+=$1} END {print sum/NR}')

    # Median calculation
    MIDDLE_INDEX=$((${#MAP_QUALS[@]}/2))
    if [ $((${#MAP_QUALS[@]} % 2)) -eq 0 ]; then
        MEDIAN=$(echo "scale=2; (${MAP_QUALS[$MIDDLE_INDEX]} + ${MAP_QUALS[$MIDDLE_INDEX-1]}) / 2" | bc)
    else
        MEDIAN=${MAP_QUALS[$MIDDLE_INDEX]}
    fi

    # Range calculation
    MIN=${MAP_QUALS[0]}
    MAX=${MAP_QUALS[@]: -1}
    RANGE=$((MAX - MIN))

    # Print the results
    echo "Mapping Quality Metrics:"
    echo "Mean: $MEAN"
    echo "Median: $MEDIAN"
    echo "Range: $RANGE (Min: $MIN, Max: $MAX)"

    # Percentage of reads above certain MAPQ thresholds
    THRESHOLDS=(10 20 30 40 50 60)
    TOTAL_READS=${#MAP_QUALS[@]}

    echo "Percentage of reads with MQ >= threshold:"
    for MAPQ in "${THRESHOLDS[@]}"; do
        ABOVE_THRESHOLD=$(printf "%s\n" "${MAP_QUALS[@]}" | awk -v threshold="$MAPQ" '$1 >= threshold' | wc -l | tr -d ' ')
        PERCENTAGE=$(awk -v above="$ABOVE_THRESHOLD" -v total="$TOTAL_READS" 'BEGIN {printf "%.2f", (above/total)*100}')
        echo "MAPQ $MAPQ: $PERCENTAGE% ($ABOVE_THRESHOLD reads)"
    done

echo ""
echo_blue "Compiling anellovirus presence data for $FULL_SAMPLE_NAME"
process_species_presence "$INPUT_BAM" "${SAMPLE_DIR}/${FULL_SAMPLE_NAME}_anellovirus_presence.txt" "${SAMPLE_DIR}/refs_with_min_5_reads.txt"

# # Generate consensus sequences for species identified
# echo_blue "Generating consensus sequences for identified species in $FULL_SAMPLE_NAME"
# samtools faidx "$ANELLOVIRUS_REFERENCEFAI"
# generate_consensus "$INPUT_BAM" "${SAMPLE_DIR}/refs_with_min_5_reads.txt" "$SAMPLE_DIR" "$FULL_SAMPLE_NAME" "$ANELLOVIRUS_REFERENCE"

# Clean up intermediate files
 echo_blue "Cleaning up intermediate files for $FULL_SAMPLE_NAME"
rm "$ANELLOMAP_FILTERED"
rm "$ANELLOMAP_OUTPUT"
rm "${SAMPLE_DIR}/${FULL_SAMPLE_NAME}_anellovirus_short.bam"
rm "${SAMPLE_DIR}/${FULL_SAMPLE_NAME}_anellovirus_long.bam"
rm "${SAMPLE_DIR}/${FULL_SAMPLE_NAME}_all_passing_reads.short.fastq.gz"
rm "${SAMPLE_DIR}/${FULL_SAMPLE_NAME}_all_passing_reads.long.fastq.gz"

# Pipeline complete
done
echo ""
echo_green "All samples have been successfully processed."

# Parse and compile the final output tables
echo ""
echo_blue "Compiling final output tables for all samples..."

# Define the path to the parsing script (same directory as this script)
PARSER_SCRIPT_PATH="$(dirname "$0")/parse_anelloflow_output.sh"

if [ -f "$PARSER_SCRIPT_PATH" ]; then
    bash "$PARSER_SCRIPT_PATH" "$INPUT_DIRECTORY"
    echo_green "Final output tables have been compiled and printed to the terminal."
else
    echo_red "Error: Parsing script not found in the same directory as Anelloflow.sh"
    echo_red "Ensure parse_anelloflow_output.sh is located alongside Anelloflow.sh."
fi
