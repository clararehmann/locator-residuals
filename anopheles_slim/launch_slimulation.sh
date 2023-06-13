#!/bin/bash

for sigma in 20 25 50; do
    mkdir -p out/slim/sigma_$sigma
    for N in {0..9}; do
        sbatch run_slimulation.batch out/slim/sigma_$sigma\/$N $sigma
    done
done
