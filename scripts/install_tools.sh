#!/bin/bash
# Descripción: Script para instalar todas las herraimentas necesarias durante el proyecto
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 27/11/2023
# Versión : 2.0

tools_dir="$HOME/tools"
mkdir -p $tools_dir
cd $tools_dir

# Function to check if a command is available
command_exists() {
  command -v "$1" >/dev/null 2>&1
}


# Verificar si el SRA Toolkit está instalado
if ! command_exists "prefetch"; then
    echo "SRA Toolkit no está instalado. Iniciando la instalación en $tools_dir..." 
    # Descargar e instalar SRA Toolkit
    echo pwd
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
    tar -xzf sratoolkit.current-centos_linux64.tar.gz
    export PATH=$PATH:$tools_dir/sratoolkit.current-centos_linux64/bin
    
    echo "SRA Toolkit instalado correctamente."  

fi

# Verificar si FastQC está instalado
if ! command_exists "fastqc" ; then
    echo "FastQC no está instalado. Iniciando la instalación en $tools_dir..."
    
    # Descargar e instalar FastQC
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
    unzip fastqc_v0.11.9.zip
    chmod +x FastQC/fastqc
    export PATH=$PATH:$tools_dir/FastQC
    
    echo "FastQC instalado correctamente."

fi


# Verificar si conda está instalado
if ! command_exists "conda"; then
    echo "Installing Miniconda..."
    mkdir -p $tools_dir/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $tools_dir/miniconda3/miniconda.sh
    bash $tools_dir/miniconda3/miniconda.sh -b -u -p $tools_dir/miniconda3
    rm -rf $tools_dir/miniconda3/miniconda.sh

    echo "Initializing conda for bash and zsh..."
    $tools_dir/miniconda3/bin/conda init bash
    $tools_dir/miniconda3/bin/conda init zsh

    echo "Updating PATH..."
    export PATH=$PATH:$tools_dir/miniconda3/condabin/
    export PATH=$PATH:$tools_dir/miniconda3/pkgs/bwa-meme-1.0.6-hdcf5f25_2/bin/

    echo "Miniconda installed and configured successfully."
else
    echo "Miniconda is already installed. Skipping installation."
fi

# Function to install a package with conda if not already installed
install_with_conda() {
  if ! command_exists "$1"; then
    echo "Installing $1 with conda..."
    conda install "$1"
  else
    echo "$1 is already installed. Skipping installation."
  fi
}

# Check and install dependencies
install_with_conda  "bedtools"
install_with_conda "samtools"
install_with_conda "bwa-meme"
install_with_conda "dragmap"

if ! command_exists "git"; then
    sudo apt install git
fi

# Function to downland git reposirtoty if not already installed
downland_with_git() {
  if ! command_exists "$1"; then
    echo "Installing $1 with ..."
    git clone $2 $tools_dir/$1
  else
    echo "$1 is already installed. Skipping installation."
  fi
}

downland_with_git "minimap2" "https://github.com/lh3/minimap2"



