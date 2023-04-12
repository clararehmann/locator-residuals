#!/bin/bash

pf7path=/sietch_colab/crehmann/pf7
#mkdir -p $pf7path
#wget -P $pf7path https://www.malariagen.net/sites/default/files/Pf7_samples.txt
#wget -P $pf7path ftp://ngs.sanger.ac.uk/production/malaria/Resource/34/Pf7.zarr.zip
unzip $pf7path/Pf7.zarr.zip -d $pf7path/genotypes_zarr
