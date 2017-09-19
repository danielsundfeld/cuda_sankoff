#!/bin/bash
BIN="../bin/cuda_sankoff"

for file in ../seqs/RFAM/*; do
    /usr/bin/time -v $BIN $file > ${file}.gpu.output 2>&1
done
