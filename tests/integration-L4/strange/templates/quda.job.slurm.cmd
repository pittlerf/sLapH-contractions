#!/bin/bash -x
#SBATCH --job-name=JOBNAME
#SBATCH --mail-type=ALL
#SBATCH --mail-user=werner@hiskp.uni-bonn.de
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:4
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --partition=batch

rundir=RUNDIR
exe=EXEC
outfile=OUTFILE
infile=INFILE
export QUDA_RESOURCE_PATH=QUDA_RSC_PATH

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

