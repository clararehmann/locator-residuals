#!/bin/bash

pf7path=/sietch_colab/crehmann/pf7
#mkdir -p $pf7path
#wget -P $pf7path https://www.malariagen.net/sites/default/files/Pf7_samples.txt
for chrom in 04 05 06 07 08 09 10 11 12 13 14; do
    wget -P $pf7path ftp://ngs.sanger.ac.uk/production/malaria/Resource/34/Pf7_vcf/Pf3D7_$chrom\_v3.pf7.vcf.gz &
done
wait
