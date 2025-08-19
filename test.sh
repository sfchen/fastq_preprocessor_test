#!/bin/bash

# Configuration
SAMPLES=("S1" "S2" "S3" "S4" "S5")  # Sample names
ADAPTER_FILE="TruSeq3-PE.fa"
R1_ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  # Illumina R1 adapter
R2_ADAPTER="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"  # Illumina R2 adapter
THREADS=4  # Adjust based on your CPU cores

# Create output directories
mkdir -p fastp_out trimmomatic_out cutadapt_out fastqc_reports logs

# Function to log timing
log_time() {
    local tool=$1
    local sample=$2
    local start=$3
    local end=$4
    echo "$(date): $tool processed $sample in $((end - start)) seconds" >> logs/timing.log
}

# Process each sample
for sample in "${SAMPLES[@]}"; do
    # Input files
    r1_input="${sample}_R1.fastq.gz"
    r2_input="${sample}_R2.fastq.gz"
    
    # Skip if input files missing
    if [ ! -f "$r1_input" ] || [ ! -f "$r2_input" ]; then
        echo "Warning: Missing files for $sample. Skipping..."
        continue
    fi

    #######################################
    # 1. Process with fastp
    #######################################
    echo "Processing $sample with fastp..."
    start=$(date +%s)
    
    fastp \
        -i "$r1_input" \
        -I "$r2_input" \
        -o "fastp_out/${sample}_R1_fastp.fastq.gz" \
        -O "fastp_out/${sample}_R2_fastp.fastq.gz" \
        -h "fastp_out/${sample}_report.html" \
        -j "fastp_out/${sample}_report.json" \
        --thread "$THREADS"
    
    end=$(date +%s)
    log_time "fastp" "$sample" "$start" "$end"

    #######################################
    # 2. Process with Trimmomatic
    #######################################
    echo "Processing $sample with Trimmomatic..."
    start=$(date +%s)
    
    trimmomatic PE \
        -threads "$THREADS" \
        -phred33 \
        "$r1_input" \
        "$r2_input" \
        "trimmomatic_out/${sample}_R1_trimmomatic_paired.fastq.gz" \
        "trimmomatic_out/${sample}_R1_trimmomatic_unpaired.fastq.gz" \
        "trimmomatic_out/${sample}_R2_trimmomatic_paired.fastq.gz" \
        "trimmomatic_out/${sample}_R2_trimmomatic_unpaired.fastq.gz" \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    end=$(date +%s)
    log_time "Trimmomatic" "$sample" "$start" "$end"

    #######################################
    # 3. Process with cutadapt
    #######################################
    echo "Processing $sample with cutadapt..."
    start=$(date +%s)
    
    cutadapt \
        -j "$THREADS" \
        -a "$R1_ADAPTER" \
        -A "$R2_ADAPTER" \
        --nextseq-trim=20 \
        -o "cutadapt_out/${sample}_R1_cutadapt.fastq.gz" \
        -p "cutadapt_out/${sample}_R2_cutadapt.fastq.gz" \
        "$r1_input" \
        "$r2_input"
    
    end=$(date +%s)
    log_time "cutadapt" "$sample" "$start" "$end"

done

#######################################
# Run FastQC on all processed files
#######################################
echo "Running FastQC on processed files..."
find  fastp_out/ trimmomatic_out/ cutadapt_out/ -name "*.fastq.gz" -print0 | \
    xargs -0 -n 1 -P "$THREADS" fastqc -o fastqc_reports/

echo "All processing complete! Results in: fastp_out/, trimmomatic_out/, cutadapt_out/"
echo "FastQC reports in: fastqc_reports/"
echo "Timing data in: logs/timing.log"
