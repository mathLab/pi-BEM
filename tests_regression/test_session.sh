#!/bin/bash

export PBEM_OPENMP=
reps=5
cooldown=30
export solvers="direct fma"

comment="
#explicit base: can only do 1-component
tag=base
echo '====== checkout ${tag} ======'
git checkout ${tag} ../include
git checkout ${tag} ../source
    
rm -rf logs_${tag}
mkdir -p logs_${tag}

#only use single-component prms
export do_reference=
export do_simple=
export do_complex=1


rm reg* ref*.log
for ((i=0; i<reps; ++i))
do 
    ./test.sh
done

mv *_log.log logs_${tag}/
"

#==============================
#these guys instead can do all
export do_reference=
export do_simple=
export do_complex=1

for tag in complex_vectors #opt001 opt002
do
    echo "====== checkout ${tag} ======"
    git checkout ${tag} ../include
    git checkout ${tag} ../source

    rm -rf logs_${tag}
    mkdir -p logs_${tag}

    rm reg* ref*.log
    for ((i=0; i<reps; ++i))
    do 
        ./test.sh
    done

    mv *_log.log logs_${tag}/
done

#explicit opt003
tag=opt003
echo "====== checkout ${tag} ======"
git checkout ${tag} ../include
git checkout ${tag} ../source
    
for omp in 0 1
do
    if (( omp != 0 ))
    then
        export PBEM_OPENMP=1
        thread="omp"
    else
        export PBEM_OPENMP=
        thread="tbb"
    fi
    
    rm -rf logs_${tag}_${thread}
    mkdir -p logs_${tag}_${thread}

    rm reg* ref*.log
    for ((i=0; i<reps; ++i))
    do 
        ./test.sh
    done

    mv *_log.log logs_${tag}_${thread}/
done

git checkout HEAD ../include
git checkout HEAD ../source

