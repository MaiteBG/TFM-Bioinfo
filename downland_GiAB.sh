#!/bin/bash

# Nombre del proyecto
PROJECT_NAME="SRR2052337"

# Directorio de descarga para el proyecto
PROJECT_DIR=$HOME/TFM/GIAB/$PROJECT_NAME

# Crear el directorio de descarga si no existe
mkdir -p $PROJECT_DIR

# Descargar el conjunto de datos del proyecto
prefetch --max-size 100G -O $PROJECT_DIR $PROJECT_NAME

# Convertir a formato FASTQ (opcional)
fastq-dump --split-files -O $PROJECT_DIR $PROJECT_NAME

# Borrar el archivo original después de la conversión
rm $PROJECT_DIR/$PROJECT_NAME.sra

# Generar informes de calidad con FastQC
fastqc -o $PROJECT_DIR $PROJECT_DIR/*.fastq


