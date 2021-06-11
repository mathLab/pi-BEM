#!/bin/bash

spack env activate wass

startdir=$(pwd)

cd ..
mkdir -p build
cd build

cmake ..
cmake --build . -j4

cp ${startdir}/ref_parameters_bem_3.prm ./parameters_bem_3.prm

for np in 1 #2 1
do 
	perf stat --detailed mpirun --np ${np} bem_fma_3d > log 2>&1
#| grep -i "we have a tria of"

    for fn in result_scalar_results scalar_error result_vector_results vector_error
    do
        cp ./${fn}.vtu ${startdir}/reg_np${np}_${fn}.vtu
    done
    
    #awk, grep or cp from log to the startdir for the info we're interested in
    cp log ${startdir}/reg_np${np}_log.log
done

cp ${startdir}/ref_parameters_bem_3_complex.prm ./parameters_bem_3.prm

for np in 1 #2 1
do 
	perf stat --detailed mpirun --np ${np} bem_fma_3d > log 2>&1
#| grep -i "we have a tria of"

    #for fn in result_scalar_results scalar_error result_vector_results vector_error
    #do
    #    cp ./${fn}.vtu ${startdir}/reg_np${np}_complex_${fn}.vtu
    #done
    
    #awk, grep or cp from log to the startdir for the info we're interested in
    cp log ${startdir}/reg_np${np}_complex_log.log
done
