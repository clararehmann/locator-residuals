#!/bin/bash

snakemake --unlock
snakemake --profile slurm --conda-frontend conda --use-conda $1
