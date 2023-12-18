#!/bin/bash

# Definir variables
GATK_VERSION="4.2.3.0"
DOCKER_IMAGE="broadinstitute/gatk:${GATK_VERSION}"
DATA_DIR="/home/maitebg/Escritorio/TFM"


# Paso 1: Preprocesamiento de datos

# Antes de analizar variantes, es necesario realizar un preprocesamiento de los datos para mejorar la calidad de las llamadas de variantes. Esto implica eliminar duplicados, realizar corrección de calidad de bases y recalibrar los valores de calidad de las variantes.

# Eliminar duplicados
sudo docker run -it --rm \
  -v "${DATA_DIR}:/datos" \
  -w /datos \
  "${DOCKER_IMAGE}" \
  gatk MarkDuplicates -I results/BWA-MEME/output.bam -O marked_duplicates.bam -M marked_dup_metrics.txt



# Corrección de calidad de bases
sudo docker run -it --rm \
  -v "${DATA_DIR}:/datos" \
  -w /datos \
  "${DOCKER_IMAGE}" \
  gatk BaseRecalibrator -I marked_duplicates.bam -R T2T-CHM13v2.0_genome/chr21.fasta -O recal_data.table

# Aplicar recalibración
sudo docker run -it --rm \
  -v "${DATA_DIR}:/datos" \
  -w /datos \
  "${DOCKER_IMAGE}" \
  gatk ApplyBQSR -I marked_duplicates.bam -R reference.fasta --bqsr-recal-file recal_data.table -O recalibrated.bam


#Paso 2: Llamada de variantes
#Una vez que los datos han sido preprocesados, puedes realizar la llamada de variantes.

# Llamada de variantes
sudo docker run -it --rm \
  -v "${DATA_DIR}:/datos" \
  -w /datos \
  "${DOCKER_IMAGE}" \
  gatk BedToIntervalList -I target.bed -O target.interval_list -SD reference.dict

sudo docker run -it --rm \
  -v "${DATA_DIR}:/datos" \
  -w /datos \
  "${DOCKER_IMAGE}" \
  gatk HaplotypeCaller -R reference.fasta -I recalibrated.bam -L target.interval_list -O raw_variants.vcf



