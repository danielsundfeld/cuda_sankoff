#!/bin/sh
nvidia-smi -a > rf03121_03124_1.fasta.out 2>&1
/usr/bin/time -v ./bin/cuda_sankoff ./seqs/sars_cov_2/rf03121_03124/rf03121_03124_1.fasta >> rf03121_03124_1.fasta.out 2>&1
aws s3 cp rf03121_03124_1.fasta.out  s3://cudasankoffresults/

/usr/bin/time -v ./bin/cuda_sankoff ./seqs/sars_cov_2/RF03121/20.fasta >> RF03121_20.fasta.out 2>&1
aws s3 cp RF03121_20.fasta.out s3://cudasankoffresults/

/usr/bin/time -v ./bin/cuda_sankoff seqs/sars_cov_2/RF03124/30.fasta >> RF03124_30.fasta.out 2>&1
aws s3 cp RF03124_30.fasta.out s3://cudasankoffresults/
sync
sudo shutdown -h now
