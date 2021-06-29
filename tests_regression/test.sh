#!/bin/bash

spack env activate wass

startdir=$(pwd)
cooldown=60

source test_utils_rc

cd ..
mkdir -p build
cd build

cmake ..
cmake --build . -j4

result=${?}

if (( result != 0 ))
then
    echo "error during build"
    exit 1
fi

#experiments (to be done both in direct solve and fma; varying mpi processes)
#-- plain scalar phi, v1
#-- plain scalar phi, v2; different function
#-- basic 2-component phi, (v1, v2)

nprocs=1
nthreads=1
ntotalcores=4

#solvers="direct fma"
solvers="fma"

do_reference=
do_simple=
do_complex=1

if [ -n "$do_reference" ]
then
    echo "reference problems (nprocs=1, nthreads=1)"
    #generate the baseline - this should be skipped if ref* files are already present
    #logs are saved as ref_* while all other isntances are reg_*
    for solver in ${solvers}
    do
        for func in 1 2
        do
            cp ${startdir}/ref_parameters_bem_3_f${func}_${solver}.prm ./parameters_bem_3.prm

            #references should be for np=1; np>1 are regressions
            nprocs=1
            nthreads=1
            
            echo "reference problem ${func} - nprocs ${nprocs} - nthreads ${nthreads} - solver ${solver}"
            if [ ! -f "${startdir}/ref_f${func}_${solver}.log" ]
            then
                run_simple
                
                sleep ${cooldown}
            else
                echo "already processed"
            fi
        done
    done
fi

if [ -n "$do_simple" ]
then
    #solve the simple problem instances with multiple mpi procs
    #logs are saved as reg_*
    for solver in ${solvers}
    do
        for func in 1 2
        do
            cp ${startdir}/ref_parameters_bem_3_f${func}_${solver}.prm ./parameters_bem_3.prm

            #references should be for np=1; np>1 are regressions
            for nprocs in 2 4
            do 
                for ((nthreads=1; nthreads <= ntotalcores/nprocs; nthreads *= 2))
                do
                    echo "regression problem ${func} - nprocs ${nprocs} - nthreads ${nthreads} - solver ${solver}"
                    run_simple
                
                    echo
                    sleep ${cooldown}
                done
                
                if (( nthreads / (ntotalcores/nprocs) != 2 ))
                then
                    #assume all procs equal
                    (( nthreads = ncores/nprocs ))
                    echo "regression problem ${func} - nprocs ${nprocs} - nthreads ${nthreads} - solver ${solver}"
                    run_simple
            
                    echo
                    sleep ${cooldown}
                fi
            done
        done
    done
fi


if [ -n "$do_complex" ]
then
    #solve complex problems - each instance shall generate two sets of result files, with f1 and f2 suffix to be compared with references
    for solver in ${solvers}
    do
        cp ${startdir}/ref_parameters_bem_3_complex_${solver}.prm ./parameters_bem_3.prm

        for nprocs in 1 2 4
        do 
            for ((nthreads=1; nthreads <= ntotalcores/nprocs; nthreads *= 2))
            do
                echo "complex problem - nprocs ${nprocs} - nthreads ${nthreads} - solver ${solver}"
                run_complex
            
                echo
                sleep ${cooldown}
            done
            
            if (( nthreads / (ntotalcores/nprocs) != 2 ))
            then
                #assume all procs equal
                (( nthreads = ncores/nprocs ))
                echo "complex problem - nprocs ${nprocs} - nthreads ${nthreads} - solver ${solver}"
                run_complex
            
                echo
                sleep ${cooldown}
            fi
        done
    done
fi

cd ${startdir}
./test_compares.sh
