#! /bin/bash
# Descripción: Ejecuta minimap2 y recoger métricas de tiempo y memoria
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 26/12/2023
# Versión : 3.0


TFM_DIR=./
RESULTS_DIR=${TFM_DIR}results/minimap/

LOG_FILE=${RESULTS_DIR}log.txt
METRICS_FILE=${RESULTS_DIR}metrics.txt

GENOME_FILE=${TFM_DIR}T2T-CHM13v2.0_genome/chr21.fasta
FASTQ_DIR=${TFM_DIR}/GIAB/

R1_FASTQ=${FASTQ_DIR}SRR2052337_1.fastq
R2_FASTQ=${FASTQ_DIR}SRR2052337_2.fastq

# Crear el directorio de resultados
mkdir -p $RESULTS_DIR

# Etiqueta para el alineamiento
echo "### Alineamiento con Minimap2" >> ${LOG_FILE}
# Realizar el alineamiento con Minimap2
{ TIME_RESULT=$( /usr/bin/time -f "%e %M" minimap2 -ax sr $GENOME_FILE $R1_FASTQ $R2_FASTQ > ${RESULTS_DIR}output.sam 2>&1 ); } 2>&1
echo "$TIME_RESULT" >> ${LOG_FILE}

${TFM_DIR}/scripts/getSummary.sh $LOG_FILE $METRICS_FILE
