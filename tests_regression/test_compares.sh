#!/bin/bash

atol=0.00001
echo "regressions instances"
for mode in direct fma
do
    for func in 1 2
    do
        for nprocs in 2 4
        do
            for nthreads in 1 2 4
            do
                for type in result_scalar_results result_vector_results #scalar_error vector_error
                do
                    if [ -f reg_np${nprocs}_nt${nthreads}_${type}_f${func}_${mode}.vtu ]
                    then
                        python mesh_compare.py ref_${type}_f${func}_${mode}.vtu reg_np${nprocs}_nt${nthreads}_${type}_f${func}_${mode}.vtu ${atol}
                        echo
                    fi
                done
            done
        done
    done
done

echo "complex instances"
for mode in direct fma
do
    for func in 1 2
    do
        for nprocs in 1 2 4
        do
            for nthreads in 1 2 4
            do
                for type in result_scalar_results #result_vector_results #scalar_error vector_error
                do
                    if [ -f reg_np${nprocs}_nt${nthreads}_${type}_complex_f${func}_${mode}.vtu ]
                    then
                        python mesh_compare.py ref_${type}_f${func}_${mode}.vtu reg_np${nprocs}_nt${nthreads}_${type}_complex_f${func}_${mode}.vtu ${atol}
                        echo
                    fi
                done
            done
        done
    done
done
