#! /bin/bash
# Descripción: Ejecuta minimap
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 27/11/2023
# Versión : 1.0



# Define variables
TFM_DIR=./
RESULTS_DIR=${TFM_DIR}results/minimap/

GENOME_FILE=${TFM_DIR}T2T-CHM13v2.0_genome/chr21.fasta
FASTQ_DIR=${TFM_DIR}/GIAB/

R1_FASTQ=${FASTQ_DIR}SRR2052337_1.fastq
R2_FASTQ=${FASTQ_DIR}SRR2052337_2.fastq

mkdir $RESULTS_DIR
minimap2 -ax sr $GENOME_FILE $R1_FASTQ $R2_FASTQ > ${RESULTS_DIR}output.sam     # paired-end alignment
