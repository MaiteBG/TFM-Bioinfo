#! /bin/bash

# Descripción: Script para seleccionar la sequencia de un cromosoma determinado dentro de un fichero fasta con todos los cromosomas (es este caso el 21)
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 23/11/2023
# Versión : 3.0

# Directorio de descarga para la sequencia fasta 
GENOME_DIR="$HOME/Escritorio/TFM/T2T-CHM13v2.0_genome"

mkdir $GENOME_DIR

#Homo sapiens isolate CHM13 chromosome 21, alternate assembly T2T-CHM13v2.0
# NCBI Reference Sequence: NC_060945.1
wget -O $GENOME_DIR/chr21.fasta  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_060945.1&rettype=fasta&retmode=text"

# Modificamos cabecera FASTA
sed -i '1s/.*/>chr21/' $GENOME_DIR/chr21.fasta 


# create index file (fai)
samtools faidx $GENOME_DIR/chr21.fasta 


