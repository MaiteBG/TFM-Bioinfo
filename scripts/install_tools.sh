#!/bin/bash
# Descripción: Script para instalar todas las herraimentas necesarias durante el proyecto
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 27/11/2023
# Versión : 2.0

tools_dir="$HOME/tools"
mkdir -p $tools_dir
cd $tools_dir

# Functio para comprobar si comando existe
command_exists() {
  command -v "$1" >/dev/null 2>&1
}


#### Conda como medio para install y compilar: bedtools, samtools, bwa-meme y dragmap ####
# Verificar si conda está instalado
if ! command_exists "conda"; then
    echo "Installing Miniconda..."
    mkdir -p $tools_dir/miniconda3
    wget https://repo.anaconda. com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $tools_dir/miniconda3/miniconda.sh
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

# Function to install a package from bioconda chanel with conda if not already installed
install_with_conda() {
  if ! command_exists "$1"; then
    echo "Installing $1 with conda..."
    conda install -c bioconda "$1"
  else
    echo "$1 is already installed. Skipping installation."
  fi
}

# Instalar paquetes a través de conda
install_with_conda "sra-tools"
install_with_conda "fastqc"
install_with_conda  "bedtools"
install_with_conda "samtools"
install_with_conda "bwa-meme"
install_with_conda "dragmap"
install_with_conda "minimap2"


# Install varaint callers

sudo apt -y update
sudo apt-get -y install docker.io

# Deepvaraint
DEEPVAR_VERSION="1.6.0"
sudo docker pull google/deepvariant:"${DEEPVAR_VERSION}"

# Install gatk
#GATK
GATK_VERSION="4.2.3.0"
sudo docker pull broadinstitute/gatk:"${GATK_VERSION}"






