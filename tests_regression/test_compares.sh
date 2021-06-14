#!/bin/bash

atol=0.005
echo "regressions instances"
for mode in direct fma
do
    for func in 1 2
    do
        for np in 2 4
        do
            for type in result_scalar_results result_vector_results #scalar_error vector_error
            do
                python mesh_compare.py ref_${type}_f${func}_${mode}.vtu reg_np${np}_${type}_f${func}_${mode}.vtu ${atol}
                echo
            done
        done
    done
done

echo "complex instances"
for mode in direct fma
do
    for func in 1 2
    do
        for np in 1 2 4
        do
            for type in result_scalar_results #result_vector_results #scalar_error vector_error
            do
                python mesh_compare.py ref_${type}_f${func}_${mode}.vtu reg_np${np}_${type}_complex_f${func}_${mode}.vtu ${atol}
                echo
            done
        done
    done
done
