#!/bin/bash

for N in 0 1 2 3 4 5 6 7 8; do
    for BW in 50000; do
        for L in 100 500; do
            zarr=out/sigma_13/run_$N\/$N\_0.05scale.zarr
            outpath=out/sigma_13/run_$N\/locator_weighted/bandwidth_$BW\_lambda_$L
            if [ -d $outpath ]; then :
            else
                mkdir -p $outpath
                for ts in {0..9}; do
                    training_set=out/sigma_13/run_$N\/training_sets/$N\_0.05scale_training_set_$ts\.txt
                    sbatch run_locator_weighted.batch $training_set $zarr $outpath\/ts_$ts $L $BW
                done
            fi
        done
    done
done

