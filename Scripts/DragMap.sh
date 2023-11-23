#! /bin/bash
# Descripción: Script para descagar y configura el entono de ejecución en caso de no estar configurado y ejecuta la busqueda para DragMap
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 16/11/2023
# Versión : 2.0


# Install and config conda if do not exist 
if ! command -v conda &> /dev/null; then
    echo "Installing Miniconda..."
    mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm -rf ~/miniconda3/miniconda.sh

    echo "Initializing conda for bash and zsh..."
    ~/miniconda3/bin/conda init bash
    ~/miniconda3/bin/conda init zsh

    echo "Updating PATH..."
    export PATH=$PATH:~/miniconda3/condabin/
    export PATH=$PATH:~/miniconda3/pkgs/bwa-meme-1.0.6-hdcf5f25_2/bin/

    echo "Miniconda installed and configured successfully."
else
    echo "Miniconda is already installed. Skipping installation."
fi

# Install bwa-meme with conda if not already installed
if ! conda list bwa-meme &> /dev/null; then
    echo "Installing bwa-meme..."
    conda install -c conda-forge -c bioconda dragmap

    echo "bwa-meme installed successfully."
else
    echo "bwa-meme is already installed. Skipping installation."
fi


# Define variables
TFM_DIR=../
GENOME_FILE=${TFM_DIR}T2T-CHM13v2.0_genome/chromosome_21.fasta
FASTQ_DIR=${TFM_DIR}GiAB/NA12878/09252015_U0a/
R1_FASTQ=${FASTQ_DIR}fastq_R1.fastq.gz
R2_FASTQ=${FASTQ_DIR}fastq_R2.fastq.gz




# Build hashtables
# We recommend using at least 8 threads
dragen-os --build-hash-table true --ht-reference $GENOME_FILE  --output-directory ${TFM_DIR}T2T-CHM13v2.0_genome/

# Run code below to train P-RMI, suffix array is required which is generated in index build code
# Takes about 15 minutes for the human genome with a single thread
build_rmis_dna.sh $GENOME_FILE



dragen-os -r ${TFM_DIR}T2T-CHM13v2.0_genome/ -1 $R1_FASTQ -2 $R2_FASTQ >  output_dragmap.sam
