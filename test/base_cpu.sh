#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH --ntasks-per-node 28
#SBATCH -t 5:00:00
SEQ="__SEQ__"
CMD="__CMD__"
OPT="__OPT__"
OUT="__OUT__"

module avail cuda
module load cuda
set -x

#run GPU program
cd $HOME"/hpc_foldalign"

strace -ve wait4 /usr/bin/time -v $CMD $OPT $SEQ >> $SEQ.$OUT.output 2>&1
