#! /bin/bash
# Descripción: Ejecutar BWA-MEME y recoger métricas de tiempo y memoria
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 26/12/2023
# Versión : 4.0

TFM_DIR=./
RESULTS_DIR=${TFM_DIR}results/BWA-MEME/

LOG_FILE=${RESULTS_DIR}log.txt
METRICS_FILE=${RESULTS_DIR}metrics.txt


GENOME_FILE=${TFM_DIR}T2T-CHM13v2.0_genome/chr21.fasta
FASTQ_DIR=${TFM_DIR}/GIAB/


THREADS=16

R1_FASTQ=${FASTQ_DIR}SRR2052337_1.fastq
R2_FASTQ=${FASTQ_DIR}SRR2052337_2.fastq

# Crear el directorio de resultados
mkdir -p $RESULTS_DIR

# Etiqueta para el índice
echo "### Construir el índice con $THREADS hilos" >> ${LOG_FILE}
# Construir el índice (toma aproximadamente 1 hora para el genoma humano)
{ TIME_RESULT=$( /usr/bin/time -f "%e %M" bwa-meme index -a meme $GENOME_FILE -t $THREADS 2>&1 ); } 2>&1
echo "$TIME_RESULT" >> ${LOG_FILE}

# Etiqueta para entrenar P-RMI
echo "### Entrenar P-RMI" >> ${LOG_FILE}
# Ejecutar código para entrenar P-RMI, se requiere un array de sufijos que se genera en el código de construcción del índice
# Toma aproximadamente 15 minutos para el genoma humano con un solo hilo
{ TIME_RESULT=$( /usr/bin/time -f "%e %M" build_rmis_dna.sh $GENOME_FILE 2>&1 ); } 2>&1
echo "$TIME_RESULT" >> ${LOG_FILE}

# Etiqueta para el alineamiento
echo "### Alineamiento con $THREADS hilos" >> ${LOG_FILE}
# Realizar el alineamiento con BWA-MEME, agregar la opción -7
{ TIME_RESULT=$( /usr/bin/time -f "%e %M" bwa-meme mem -7 -Y -K 100000000 -t $THREADS $GENOME_FILE $R1_FASTQ $R2_FASTQ -o ${RESULTS_DIR}output_${THREADS}.sam 2>&1 ); } 2>&1
echo "$TIME_RESULT" >> ${LOG_FILE}



${TFM_DIR}/scripts/getSummary.sh $LOG_FILE $METRICS_FILE
