#! /bin/bash
# Descripción: Ejecutar BWA-MEME y recoger métricas de tiempo
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 18/12/2023
# Versión : 3.0

#!/bin/bash

TFM_DIR=./
RESULTS_DIR=${TFM_DIR}results/BWA-MEME/
GENOME_FILE=${TFM_DIR}T2T-CHM13v2.0_genome/chr21.fasta
FASTQ_DIR=${TFM_DIR}/GIAB/

R1_FASTQ=${FASTQ_DIR}SRR2052337_1.fastq
R2_FASTQ=${FASTQ_DIR}SRR2052337_2.fastq

# Crear el directorio de resultados
mkdir -p $RESULTS_DIR

# Guardar la fecha y hora de inicio
start_time=$(date +"%Y-%m-%d %H:%M:%S")


# Construir el índice (toma aproximadamente 1 hora para el genoma humano)
# Recomendamos usar al menos 8 hilos

index_start_time=$(date +"%Y-%m-%d %H:%M:%S")
bwa-meme index -a meme $GENOME_FILE -t 8
index_end_time=$(date +"%Y-%m-%d %H:%M:%S")


# Ejecutar código para entrenar P-RMI, se requiere un array de sufijos que se genera en el código de construcción del índice
# Toma aproximadamente 15 minutos para el genoma humano con un solo hilo
build_start_time=$(date +"%Y-%m-%d %H:%M:%S")
build_rmis_dna.sh $GENOME_FILE
build_end_time=$(date +"%Y-%m-%d %H:%M:%S")

# Realizar el alineamiento con BWA-MEME, agregar la opción -7
alignment_start_time=$(date +"%Y-%m-%d %H:%M:%S")
bwa-meme mem -7 -Y -K 100000000 -t 8 $GENOME_FILE $R1_FASTQ $R2_FASTQ -o ${RESULTS_DIR}output.sam
alignment_end_time=$(date +"%Y-%m-%d %H:%M:%S")


# Guardar la fecha y hora de finalización
end_time=$(date +"%Y-%m-%d %H:%M:%S")

# Calcular y guardar los tiempos de ejecución en segundos

end_seconds=$(date -d "$end_time" '+%s')

echo -e "Tiempo_total\t$((end_seconds - start_seconds))" > ${RESULTS_DIR}metrics.txt
echo "Tiempo_construir_índice\t$((index_start_time - index_end_time))" >> ${RESULTS_DIR}metrics.txt
echo "Tiempo_entrenar_P-RMI\t$((build_start_time - build_end_time))" >> ${RESULTS_DIR}metrics.txt
echo "Tiempo_alineamiento\t$((alignment_start_time - alignment_end_time))" >> ${RESULTS_DIR}metrics.txt

