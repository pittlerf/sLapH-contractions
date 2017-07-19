#!/bin/bash
maxJobN=60 # maximum of jobs in queue
name=sWC_A2p1

# loop for the configs
for cfg in {501..549..8}; do
  echo "submitting conf ${cfg}"
 
    ## check whether at maximum of jobs in queue
    while [ $(qstat -u werner | grep -v "R" | grep $name | wc -l) -gt $maxJobN ]; do
      echo Waiting, queue is at maximum.
      sleep 30
    done
  
    j=`printf %d ${cfg}`
  
    cd cnfg${j}
    qsub job_script.pbs
    cd ../
  
done

#minSpace=-5 # TB of quota left
  ## check for quota space
#  while [ $(grep work /homec/hbn28/usage.quota |head -1|awk -v min=$minSpace '{if(int(($4-$3)/1000)>min){ print(1);} else {print(0)}}') -eq 0 ] ; do
#    echo -n "out of space! "
#    echo $(grep hbn287 /homec/hbn28/usage.quota |awk '(NR==3){print(int($3/1000))}') TB used
#    sleep 5m
#  done


#randomdir="rnd_vec_1"
#
#for i in {1000..2696..8}; do
##for i in 2048; do
#  echo "starting job on config $i"
#  j=`printf %04d ${i}`
#
#  cd cnfg${j}/${randomdir}
#
#  qsub job_script.pbs
#
#  cd ../../
#done
