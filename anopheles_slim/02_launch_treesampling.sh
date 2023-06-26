#!/bin/bash


# sigma 25

for N in {0..9}; do
    tree=out/slim/sigma_25/$N\.trees
    if [ -f $tree ]; then
        if [ ! -d out/sigma_25/run_$N\/training_sets ]; then
            mkdir -p out/sigma_25/run_$N\/training_sets
            sbatch run_treesampling.batch $tree out/sigma_25/run_$N\/$N\_0.05scale out/sigma_25/run_$N\/training_sets/$N\_0.05scale
        fi
    fi
done
