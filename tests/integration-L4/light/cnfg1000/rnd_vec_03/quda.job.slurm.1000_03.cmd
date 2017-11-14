#!/bin/bash -x
#SBATCH --job-name=test4x4x4x4_u_light_1000_3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=werner@hiskp.uni-bonn.de
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:4
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --partition=batch

rundir=/hiskp2/ueding/test4x4x4x4/light/cnfg1000/rnd_vec_03
exe=/hadron/bartek/bin/peram_gen/peram_gen.multigpu.hybrid.quda-v0.7.2.openmpi.ranlxd2
outfile=../outputs/run_1000_03.out
infile=LapH_1000_03.in
export QUDA_RESOURCE_PATH=/hadron/bartek/misc/quda_resources/v0.7.2-openmpi

if [ ! -d ${QUDA_RESOURCE_PATH} ]; then
  mkdir -p ${QUDA_RESOURCE_PATH}
fi

NUMACTL= #numactl -N 0 -m 0

cd ${rundir}
date > ${outfile}
QUDA_RESOURCE_PATH=${QUDA_RESOURCE_PATH} OMP_NUM_THREADS=1 \
  QUDA_ENABLE_GDR=1 QUDA_ENABLE_P2P=1 QUDA_ENABLE_TUNING=1 \
  mpirun.openmpi -np ${SLURM_NTASKS} ${exe} -LapHsin ${infile} | tee -a ${outfile}

date >> ${outfile}

