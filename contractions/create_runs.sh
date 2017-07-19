#!/bin/bash


for i in {501..549..8}; do
  echo "creating config $i"

  # create dir
  mkdir -p cnfg$i
  cd cnfg$i

  # create correct job script
  cp ../archetype/job_script.pbs .
  perl -pe "s/.*/#PBS -N sWC_A2p1_${i}  / if $. == 1" < job_script.pbs > job_script.pbs2
  perl -pe "s/.*/cd \/hiskp2\/werner\/sWC_A2p1_Mpi270_L24T96_TF1_apbc-no-compression_2\/cnfg${i}  / if $. == 8" < job_script.pbs2 > job_script.pbs
#  mv job_script.pbs2 job_script.pbs
  rm job_script.pbs2

  # adapt contract.in
  cp ../archetype/contract.in .
  perl -pe "s/.*/start_config = "$i"/ if $. == 15" < contract.in > contract.in2
  perl -pe "s/.*/end_config   = "$i"/ if $. == 16" < contract.in2 > contract.in
  rm contract.in2

#  # adapt input files
#  for j in {0..4..1}; do
#    cp ../archetype/mom$j.in .
#    perl -pe "s/.*/start_config = "$i"/ if $. == 15" < mom$j.in > mom$j.in2
#    perl -pe "s/.*/end_config   = "$i"/ if $. == 16" < mom$j.in2 > mom$j.in
#    rm mom$j.in2
#  done

  # copy executables
  cp ../archetype/contract .
  cd ../

done
