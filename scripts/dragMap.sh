#! /bin/bash
# Descripción: Ejecutar dragMap y recoger métricas de tiempo y memoria
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 26/11/2023
# Versión : 3.0



# Define variables
TFM_DIR=./
RESULTS_DIR=${TFM_DIR}results/dragMap/

LOG_FILE=${RESULTS_DIR}log.txt
METRICS_FILE=${RESULTS_DIR}metrics.txt

GENOME_FILE=${TFM_DIR}T2T-CHM13v2.0_genome/chr21.fasta
FASTQ_DIR=${TFM_DIR}GIAB/
R1_FASTQ=${FASTQ_DIR}SRR2052337_1.fastq.gz
R2_FASTQ=${FASTQ_DIR}SRR2052337_2.fastq.gz

# Crear el directorio de resultados
mkdir -p $RESULTS_DIR

# Etiqueta para la construcción de las tablas hash
echo "### Construcción de las tablas hash" >> $LOG_FILE
{ TIME_RESULT=$( /usr/bin/time -f "%e %M" dragen-os --build-hash-table true --ht-reference $GENOME_FILE --output-directory ${TFM_DIR}T2T-CHM13v2.0_genome/ 2>&1 ); } 2>&1
echo "$TIME_RESULT" >> $LOG_FILE


# Etiqueta para el alineamiento con dragmap
echo "### Alineamiento con dragmap" >> $LOG_FILE
# Realizar el alineamiento con dragmap
{ TIME_RESULT=$( /usr/bin/time -f "%e %M" dragen-os -r ${TFM_DIR}T2T-CHM13v2.0_genome/ -1 $R1_FASTQ -2 $R2_FASTQ > ${RESULTS_DIR}output_dragmap.sam 2>&1 ); } 2>&1
echo "$TIME_RESULT" >> $LOG_FILE



${TFM_DIR}/scripts/getSummary.sh $LOG_FILE $METRICS_FILE
