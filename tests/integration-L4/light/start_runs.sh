#!/bin/bash

conf_start=1501
conf_step=16
conf_end=1997

rv_start=0
rv_end=5

for i in $( seq ${conf_start} ${conf_step} ${conf_end} ); do
  j=`printf %04d $i`

  cd cnfg${j}/
  
  for rv in $( seq ${rv_start} ${rv_end} ); do
    echo "starting config $i rv $rv"
    rnd2=$( printf %02d ${rv} )
    cd rnd_vec_${rnd2}
  
    jscr=quda.job.pbs.${j}_${rnd2}.sh
 
    qsub ${jscr}

    cd ..
  done
  cd ..
done
