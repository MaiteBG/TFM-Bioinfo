#!/bin/bash
# Descripción: Script para tranformar un un fichero sam que cotnenga referneiscias al cromsoomas 21 a un fichero bam y bed
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 18/12/2023
# Versión : 3.0

# Input SAM file
input_sam="$1"

# Output BAM file (remove .sam extension and add .bam)
output_bam="${input_sam%.sam}.bam"

# Output BED file (remove .sam extension and add .bed)
output_bed="${input_sam%.sam}.bed"

# Convert SAM to BAM using samtools
echo "Converting $input_sam to $output_bam..."
samtools view -S -b "$input_sam" > "$output_bam"


# Convert sorted BAM to BED using bedtools
echo "Converting ${output_bam%.bam}.bam to $output_bed..."

#Si en neceario cambiar el nombre de la sequencia por el del cromosoma
bedtools bamtobed -i "${output_bam%.bam}.bam" |  awk '{print "chr21\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11}' > "$output_bed"

#Ordenamos bed
sort -k1,1 -k2,2n -o ./output.bed ./output.bed

#Comprimimos y obtenemos indices
bgzip output.bed
tabix -p bed output.bed.gz

echo "Conversion completed. BED file saved as $output_bed."

