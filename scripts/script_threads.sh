#!/bin/bash
#PBS -N bem_thread_scaling
#PBS -l walltime=12:00:00
#PBS -l nodes=2:ppn=20
#PBS -q regular




module purge

WORK_DIR=/home/ngiuliani/parallel-bem/build_clash_thread
. /home/mathlab/gnu.conf /home/mathlab/gnu/

cd $WORK_DIR

#for i in 1 2 4 8 16 20
#do
#  mpirun -np 1 --bind-to none bem_fma $i >> output_bem_fma_threads_$((${i})).txt
#done
mpirun -np 6 -pernode --bind-to none bem_fma 6 >> output_bem_fma_threads_36_6_6.txt 
mpirun -np 8 -pernode --bind-to none bem_fma 5 >> output_bem_fma_threads_40_8_5.txt
mpirun -np 4 -pernode --bind-to none bem_fma 10 >> output_bem_fma_threads_40_4_10.txt
mpirun -np 2 -pernode --bind-to none bem_fma 20 >> output_bem_fma_threads_40_2_20.txt
