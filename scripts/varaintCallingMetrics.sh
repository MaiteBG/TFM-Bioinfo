#!/bin/bash
# Descripción:
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 27/12/2023
# Versión : 1.0


ALGORITHM_DIR=$1

REF_PATH="/chain_and_TP/chr21_HG38.fasta"
DATA_DIR="/home/maitebg/Escritorio/TFM"


# Verificar si el número de parámetros es 1
if [ "$#" -eq 0 ]; then
    echo "Falta el paramentros nombre_algortimo"
    exit 1
fi

# Verificar que el directorio del algoritmo existe y no está vacío
echo ${DATA_DIR}/results/${ALGORITHM_DIR}
if [ ! -d "${DATA_DIR}/results/${ALGORITHM_DIR}" ] || [ -z "$(ls -A ${DATA_DIR}/results/${ALGORITHM_DIR})" ]; then
  echo "Error: El directorio del algoritmo no existe o está vacío."
  exit 2
fi


sudo docker pull pkrusche/hap.py


sudo docker run   \
  -v "${DATA_DIR}":"/data" \
  -e HGREF=/data${REF_PATH} \
  pkrusche/hap.py /opt/hap.py/bin/hap.py \
  /data/chain_and_TP/chr21_TP_norm.vcf.gz\
   /data/results/${ALGORITHM_DIR}/output_after_chain.vcf.gz \
  -f /data/chain_and_TP/chr21.bed.gz \
  -o /data/results/${ALGORITHM_DIR}/happy.output 
 





