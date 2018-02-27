#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH --ntasks-per-node 28
#SBATCH -t 5:00:00
SEQ="./seqs/RFAM_BP/RF00127.fasta"
CMD="./bin/sankoff"
OPT="-t 1"
OUT="cpu1"

module avail cuda
module load cuda
set -x

#run GPU program
cd $HOME"/hpc_foldalign"

strace -ve wait4 /usr/bin/time -v $CMD $OPT $SEQ >> $SEQ.$OUT.output 2>&1
