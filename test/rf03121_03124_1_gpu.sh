#!/bin/sh
nvidia-smi -a > rf03121_03124_1.fasta.out 2>&1
/usr/bin/time -v ./bin/cuda_sankoff ./seqs/sars_cov_2/rf03121_03124/rf03121_03124_1.fasta >> rf03121_03124_1.fasta.out 2>&1
aws s3 cp rf03121_03124_1.fasta.out  s3://cudasankoffresults/
sudo shutdown -h now
