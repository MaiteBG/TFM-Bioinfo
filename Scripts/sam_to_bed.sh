#!/bin/bash

# Function to check if a command is available
command_exists() {
  command -v "$1" >/dev/null 2>&1
}

# Function to install a package with conda if not already installed
install_with_conda() {
  if ! command_exists "$1"; then
    echo "Installing $1 with conda..."
    conda install -c bioconda "$1"
  else
    echo "$1 is already installed. Skipping installation."
  fi
}

# Check and install dependencies
install_with_conda "bedtools"
install_with_conda "samtools"

# Check if input SAM file is provided as a parameter
if [ "$#" -eq 0 ]; then
  echo "Usage: $0 <input.sam>"
  exit 1
fi

# Input SAM file
input_sam="$1"

# Output BAM file (remove .sam extension and add .bam)
output_bam="${input_sam%.sam}.bam"

# Output BED file (remove .sam extension and add .bed)
output_bed="${input_sam%.sam}.bed"

# Convert SAM to BAM using samtools
echo "Converting $input_sam to $output_bam..."
samtools view -S -b "$input_sam" > "$output_bam"

# Convert BAM to BED using bedtools
echo "Converting $output_bam to $output_bed..."
bedtools bamtobed -i "$output_bam" | awk '{print "chr21\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11}' > "$output_bed"

echo "Conversion completed. BED file saved as $output_bed."

