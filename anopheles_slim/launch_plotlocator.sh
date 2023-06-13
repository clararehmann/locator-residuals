#!/bin/bash

#!/bin/bash


for N in 0 1 2 3 4 5 6 7 8; do
    for BW in 50000; do
        for L in 100 500; do
            metadata=out/sigma_13/run_$N\/$N\_0.05scale_metadata.txt
            locator=out/sigma_13/run_$N\/locator_weighted/bandwidth_$BW\_lambda_$L
            sbatch run_plotlocator.batch $locator $metadata out/centroids/sigma_13_$N\_0.05scale_bandwidth_$BW\_lambda_$L
        done
    done
done


#mkdir -p out/centroids

#for N in {0..9}; do
#    locator=out/sigma_25/run_$N\/locator
#    metadata=out/sigma_25/run_$N\/$N\_0.05scale_metadata.txt
#    sbatch run_plotlocator.batch $locator $metadata out/centroids/sigma_25_$N\_0.05scale
#done
