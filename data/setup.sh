#!/bin/bash

# Get Java
sudo apt-get install openjdk-8-jdk

# Install nextflow (& move to path)
wget -qO- https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/bin/nextflow

# pull docker images
docker pull staphb/quast:5.0.2
docker pull staphb/prokka:1.14.5
docker pull staphb/ncbi-amrfinderplus:3.10.5
docker pull staphb/piggy:1.5
docker pull quay.io/biocontainers/panaroo:1.2.3--py_0
docker pull staphb/snippy:4.6.0
docker pull rocker/tidyverse:4.0.5
docker pull quay.io/biocontainers/scoary:1.6.16--py_2
docker pull staphb/snp-dists:0.8.2
