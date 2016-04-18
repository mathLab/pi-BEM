#!/bin/bash
#PBS -N bem_local_convergence
#PBS -l walltime=12:00:00
#PBS -l nodes=2:ppn=20
#PBS -q regular




module purge

WORK_DIR=/home/ngiuliani/parallel-bem/build_local_Q3
. /home/mathlab/candi-testing.conf

cd $WORK_DIR
rm -rf *Mak* *mak* *bem_fma* *txt *bin *vtu
cmake ../
make -j20

#for i in 1 2 4 8 16 20
#do
#  mpirun -np 1 --bind-to none bem_fma $i >> output_bem_fma_threads_$((${i})).txt
#done
cycle_old=2
old_str="set Number of cycles  = 2"

for i in 0 1 2 3 4 5 6 7 8 9 
do
  cycle=$i
  # cycle_str=$( printf '%04d' $cycle )                                                                                                                                            
  # cycle_old_str=$(printf '%04d' $cycle_old )                                                                                                                                     
  #new_str=sed "s/$cycle_old_str/$cycle_str/g" "$old_str"                                                                                                                          
  new_str=${old_str/$cycle_old/$cycle}
  echo $cycle_str
  echo $new_str
  f="parameter_bem.prm"
  sed "s/$old_str/$new_str/g" "$f" > pippo.txt
  cp pippo.txt parameter_bem.prm
  cycle_old=$cycle
  old_str=$new_str
  mpirun -np 2 --map-by ppr:1:node --bind-to none bem_fma >> output_bem_fma_$((${i}))_local_refs_FMM.txt 
done

