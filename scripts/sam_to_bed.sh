#!/bin/bash
# Descripción: Script para transformar un fichero SAM que contiene referencias al cromosoma 21 a un fichero BAM y BED
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 26/12/2023
# Versión: 2.0

# Comprueba si se proporciona un archivo SAM como argumento
if [ $# -eq 0 ]; then
    echo "Uso: $0 <input_sam>"
    exit 1
fi

# Input SAM file
input_sam="$1"

# Verifica si el archivo SAM de entrada existe
if [ ! -e "$input_sam" ]; then
    echo "Error: El archivo $input_sam no existe."
    exit 1
fi

# Output BAM file (remueve la extensión .sam y agrega .bam)
output_bam="${input_sam%.sam}.bam"

# Output BED file (remueve la extensión .sam y agrega .bed)
output_bed="${input_sam%.sam}.bed"

# Convertir SAM a BAM ordenado usando samtools
echo "Convirtiendo $input_sam a $output_bam..."
samtools sort -o "$output_bam" "$input_sam"
samtools index "$output_bam"

# Convertir BAM ordenado a BED usando bedtools
echo "Convirtiendo ${output_bam%.bam}.bam a $output_bed..."
# Si es necesario, cambiar el nombre de la secuencia por el del cromosoma (descomentar la línea siguiente)
# bedtools bamtobed -i "${output_bam%.bam}.bam" | awk '{print "chr21\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11}' > "$output_bed"
bedtools bamtobed -i "${output_bam%.bam}.bam" > "$output_bed"


# Ordenar el archivo BED (como BAM esta ordenador no es necesario)
# sort -k1,1 -k2,2n -o "$output_bed" "$output_bed"


# Comprimir y obtener índices
bgzip "$output_bed"
tabix -p bed "$output_bed.gz"

echo "Conversión completada. Archivo BED guardado como $output_bed."

