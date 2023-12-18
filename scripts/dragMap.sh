#! /bin/bash
# Descripción: Script para descagar y configura el entono de ejecución en caso de no estar configurado y ejecuta la busqueda para DragMap
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 16/11/2023
# Versión : 2.0



# Define variables

TFM_DIR=./
RESULTS_DIR=${TFM_DIR}results/dragMap/
GENOME_FILE=${TFM_DIR}T2T-CHM13v2.0_genome/chr21.fasta
FASTQ_DIR=${TFM_DIR}GIAB/
R1_FASTQ=${FASTQ_DIR}SRR2052337_1.fastq.gz
R2_FASTQ=${FASTQ_DIR}SRR2052337_2.fastq.gz



# Build hashtables
# We recommend using at least 8 threads
#dragen-os --build-hash-table true --ht-reference $GENOME_FILE  --output-directory ${TFM_DIR}T2T-CHM13v2.0_genome/

# Run code below to train P-RMI, suffix array is required which is generated in index build code
# Takes about 15 minutes for the human genome with a single thread
#build_rmis_dna.sh $GENOME_FILE

# Run dragmap
dragen-os -r ${TFM_DIR}T2T-CHM13v2.0_genome/ -1 $R1_FASTQ -2 $R2_FASTQ >  output_dragmap.sam


