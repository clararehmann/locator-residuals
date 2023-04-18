#!/bin/bash

snakemake -j 100 --use-conda --rerun-incomplete --keep-going $1
