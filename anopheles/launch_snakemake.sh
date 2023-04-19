#!/bin/bash

snakemake --unlock
snakemake -j 5 --use-conda --rerun-incomplete --keep-going $1
