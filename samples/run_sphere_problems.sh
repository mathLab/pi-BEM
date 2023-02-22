#!/bin/bash

mkdir -p logs

cad_refine=4

export OMP_THREADS_NUM=4

for problem in real complex
do
    for bconds in DN DR NR
    do
        for solver in direct fma
        do
            echo "problem ${problem} - boundaries ${bconds} - solver ${solver}"
            cp ref_params_sphere_${problem}_${bconds}_${solver}.prm ./parameters_bem_3.prm
            
            #tune the number of refinements
            set Number of cycles                                                         = 2
            sed -i "s/set Number of cycles                                                         = .*/set Number of cycles                                                         = ${cad_refine}/" parameters_bem_3.prm
            
            perf stat --detailed ./bem_fma_3d > log 2>&1
            
            mv log logs/sphere_${problem}_${bconds}_${solver}.log
            
            mkdir -p logs/sphere_${problem}_${bconds}_${solver}
            mv meshResult.inp *.vtu logs/sphere_${problem}_${bconds}_${solver}
        done
    done
done
