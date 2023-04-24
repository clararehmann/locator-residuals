#!/bin/bash

snakemake --unlock
snakemake --profile slurm --conda-frontend conda --use-conda --keep-going $1
