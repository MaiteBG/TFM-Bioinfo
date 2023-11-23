#! /bin/bash
# Descripción: Script para descagar y configura el entono de ejecución en caso de no estar configurado y ejecuta la busqueda para BWA-MEME
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
    conda install -c conda-forge -c bioconda bwa-meme

    echo "bwa-meme installed successfully."
else
    echo "bwa-meme is already installed. Skipping installation."
fi
# Print version and mode of compiled binary executable
bwa-meme version


# Define variables
TFM_DIR=./
ALGORTIM_DIR=
GENOME_FILE=${TFM_DIR}T2T-CHM13v2.0_genome/chromosome_21.fasta
FASTQ_DIR=${TFM_DIR}GiAB/NA12878/09252015_U0a/
R1_FASTQ=${FASTQ_DIR}fastq_R1.fastq
R2_FASTQ=${FASTQ_DIR}fastq_R2.fastq


# Build index (Takes ~1hr for human genome)
# We recommend using at least 8 threads
bwa-meme index -a meme $GENOME_FILE -t 8

# Run code below to train P-RMI, suffix array is required which is generated in index build code
# Takes about 15 minutes for the human genome with a single thread
build_rmis_dna.sh $GENOME_FILE

# Perform alignment with BWA-MEME, add -7 option
bwa-meme mem -7 -Y -K 100000000 -t 8 $GENOME_FILE $R1_FASTQ $R2_FASTQ -o ${TFM_DIR}output_meme.sam

