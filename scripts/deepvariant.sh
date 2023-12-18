#! /bin/bash
# Descripción: Ejecuta varaint caller con deepvariant a traves de un docker
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 10/12/2023
# Versión : 1.0


ALGORITHM_DIR=$1
BIN_VERSION="1.6.0"
REF_PATH="T2T-CHM13v2.0_genome/chr21.fasta"
DATA_DIR="/home/maitebg/Escritorio/TFM"

# Verificar que el directorio del algoritmo existe y no está vacío
echo ${DATA_DIR}/results/${ALGORITHM_DIR}
if [ ! -d "${DATA_DIR}/results/${ALGORITHM_DIR}" ] || [ -z "$(ls -A ${DATA_DIR}/results/${ALGORITHM_DIR})" ]; then
  echo "Error: El directorio del algoritmo no existe o está vacío."
  exit 1
fi

sudo docker run \
  -v "${DATA_DIR}":"/data" \
  -w /data \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref="${REF_PATH}" \
  --reads="results/${ALGORITHM_DIR}/output.bam" \
  --regions "NC_060945.1:10,000,000-10,010,000" \
  --output_vcf="results/${ALGORITHM_DIR}/output.vcf.gz" \
  --output_gvcf="results/${ALGORITHM_DIR}/output.g.vcf.gz" \
  --intermediate_results_dir="results/${ALGORITHM_DIR}/intermediate_results_dir" \
  --num_shards=1

