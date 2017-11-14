#!/bin/bash

conf_start=1000 
conf_step=1
conf_end=1000 

confs=$1
if [ -z "$1" ]; then
  confs=$( seq ${conf_start} ${conf_step} ${conf_end} )
fi

rvs=""
rv_start=$2
rv_end=$3
if [ -z "$2" -o -z "$3" ]; then
  rvs=$(seq 0 4)
else
  rvs=$(seq $rv_start $rv_end)
fi

for i in $confs; do
  j=`printf %04d $i`

  cd cnfg${j}/
  
  for rv in $rvs; do
    echo "starting config $i rv $rv"
    rnd2=$( printf %02d ${rv} )
    cd rnd_vec_${rnd2}
  
    jscr=quda.job.slurm.${j}_${rnd2}.cmd
 
    sbatch ${jscr}

    cd ..
  done
  cd ..
done
