#!/bin/sh
nvidia-smi -a > RF03122_51.fasta.out.txt 2>&1
/usr/bin/time -v ./bin/cuda_sankoff seqs/sars_cov_2/RF03122/51.fasta >> RF03122_51.fasta.out.txt 2>&1
aws s3 cp RF03122_51.fasta.out.txt  s3://cudasankoffresults/

/usr/bin/time -v ./bin/cuda_sankoff seqs/sars_cov_2/delta/reference_delta.fasta >> ref_delta_gpu.fasta.out.txt 2>&1
aws s3 cp ref_delta_gpu.fasta.out.txt s3://cudasankoffresults/

cat /proc/cpuinfo >> ref_delta_cpu.fasta.out.txt
/usr/bin/time -v ./bin/sankoff seqs/sars_cov_2/delta/reference_delta.fasta >> ref_delta_cpu.fasta.out.txt 2>&1
aws s3 cp ref_delta_cpu.fasta.out.txt s3://cudasankoffresults/
sync
sudo shutdown -h now
