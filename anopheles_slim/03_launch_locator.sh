#!/bin/bash



# sigma 25

#for N in {0..9}; do
#    tree=out/slim/sigma_20/$N\.trees
#    if [ -f $tree ]; then
#        if [ ! -d out/sigma_20/run_$N\/training_sets ]; then
#            mkdir -p out/sigma_20/run_$N\/training_sets
#            sbatch run_treesampling.batch $tree out/sigma_20/run_$N\/$N\_0.05scale out/sigma_20/run_$N\/training_sets/$N\_0.05scale
#        fi
#    fi
#done


for N in {0..9}; do
    zarr=out/sigma_25/run_$N\/$N\_0.05scale.zarr
    outpath=out/sigma_25/run_$N\/locator
    mkdir -p $outpath
    if [ -d $zarr ]; then
        for ts in {0..9}; do
            training_set=out/sigma_25/run_$N\/training_sets/$N\_0.05scale_training_set_$ts\.txt
            predlocs=$(ls $outpath\/out_$ts\_*_predlocs.txt)
        
        ## identify window to start at    

            window=0            
            for p in $predlocs; do 
                IFS='_' read -ra WINSIZE <<< $p
                IFS='-' read -ra WINSTOP <<< ${WINSIZE[3]} 
                winstop=${WINSTOP[1]}
                if [ $winstop -gt $window ]; then
                    window=$winstop
                fi
            done
            
            if [ $window -eq 0 ]; then :
            else window=$(expr $window + 1)
            fi
            
        ## start Locator run @ next window, if not finished running
            
            if [ $window -eq 100000000 ]; then :
            else sbatch run_locator.batch $training_set $zarr $outpath\/out_$ts $window
            fi
        done
    fi
done
