#!/bin/bash

snakemake -j 10 --use-conda --rerun-incomplete --keep-going $1
