#! /bin/bash
# Descripción: Script para descagar los ficheros .fastq de prueva el inidividuo de giab 
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 14/11/2023
# Versión : 1.0


# Descargamos indices de las sequencias de Illumina300X_wgs_09252015 de NA12878
url_index="https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_indexes/NA12878/"
index_file_name="sequence.index.NA12878_Illumina300X_wgs_09252015_updated"
directory_name="../GiAB/NA12878_09252015_U0a/"

mkdir -p $directory_name
cd $directory_name


wget $url_index$index_file_name  # Obtenemos indices


# Recoremos el fichero y descargamos los ficheros de la Sample_U0a de la update del 131219_D00360_005_BH814YADXX/Project_RM8398/Sample_U0a


gawk '$1 ~ "/131219_D00360_005_BH814YADXX/Project_RM8398/Sample_U0a/U0a_.*001.fastq.gz" { system(" wget " $1 ) ; system(" wget " $3 ) };'  $index_file_name



# Juntamos ficheros pair-end R1 y R2 (los apidend teiene el mismo nnombre solo cambiando R1/2 y zcat los ordena de forma alfabetica por lo tanto se mantiene el orden)




