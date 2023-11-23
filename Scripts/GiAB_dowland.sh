#! /bin/bash
# Descripción: Script para descagar los ficheros .fastq de prueva el inidividuo de giab 
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 14/11/2023
# Versión : 1.0


# Descargamos indices de las sequencias de Illumina300X_wgs_09252015 de NA12878
url_index="https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_indexes/NA12878/"
index_file_name="sequence.index.NA12878_Illumina300X_wgs_09252015_updated"
directory_name="./GiAB/NA12878/"

mkdir -p $directory_name
cd $directory_name


wget $url_index$index_file_name # Obtenemos indices


# Recoremos el fichero y descargamos los ficheros de la Sample_U0a del 131219

gawk '$1 ~ "/131219.*U0a_.*001.fastq.gz" { print $1 };' ./$index_file_name

# Juntamos ficheros pair-end R1 y R2

zcat ./*R1_001.fastq.gz > reads_R1.fastq

zcat ./*R2_001.fastq.gz > reads_R2.fastq



