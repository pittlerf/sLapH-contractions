#!/bin/bash

conf_step=1
conf_start=1000
conf_end=1000

rv_start=0
rv_end=4

seeds=(83456 87154 12589 73155 13426)

flavour="strange"

RUNDIR="/hiskp2/ueding/test4x4x4x4/strange/cnfg"
EVDIR="/hiskp2/eigensystems/test4x4x4x4/nev_16"
GCONFBASE="/hiskp2/gauges/test4x4x4x4/conf"
EXEC="/hadron/bartek/bin/peram_gen/peram_gen.multigpu.hybrid.quda-v0.7.2.openmpi.ranlxd2"
JOBNAME="test4x4x4x4_u_${flavour}"
QUDA_RSC_PATH="/hadron/bartek/misc/quda_resources/v0.7.2-openmpi"

for i in $( seq ${conf_start} ${conf_step} ${conf_end} ); do
  echo "creating config $i"
  j=`printf %04d $i`

  mkdir -p cnfg${j}/
  mkdir -p cnfg${j}/outputs
  cd cnfg${j}/
  
  for rv in $( seq ${rv_start} ${rv_end} ); do
    seed=${seeds[${rv}]}
    rnd2=$( printf %02d ${rv} )
    mkdir -p rnd_vec_${rnd2}
    cd rnd_vec_${rnd2}
  
    cp ../../templates/quda.invert.input invert.input
    sed -i "s@NSTORE@${i}@g" invert.input
    sed -i "s@GCONFBASE@${GCONFBASE}@g" invert.input
    
    laphin=LapH_${j}_${rnd2}.in
    jscr=quda.job.slurm.${j}_${rnd2}.cmd
    outfile="../outputs/run_${j}_${rnd2}.out"
  
    cp ../../templates/quda.job.slurm.cmd ${jscr}
    sed -i "s@RUNDIR@${RUNDIR}${j}/rnd_vec_${rnd2}@g" ${jscr}
    sed -i "s@JOBNAME@${JOBNAME}_${j}_${rv}@g" ${jscr}
    sed -i "s@INFILE@${laphin}@g" ${jscr}
    sed -i "s@OUTFILE@${outfile}@g" ${jscr}
    sed -i "s@EXEC@${EXEC}@g" ${jscr}
    sed -i "s@QUDA_RSC_PATH@${QUDA_RSC_PATH}@g" ${jscr}
  
    cp ../../templates/quda.LapH.in ${laphin}
    sed -i "s@NSTORE@${i}@g" ${laphin} 
    sed -i "s@NB_RND@${rv}@g" ${laphin}
    sed -i "s@SEED@${seed}@g" ${laphin}
    sed -i "s@EVDIR@${EVDIR}@g" ${laphin}

    cd ..
  done
  cd ..
done
