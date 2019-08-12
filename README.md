# CUDA-Sankoff

CUDA-Sankoff is a tool that performs the optimal RNA structural alignment. We developed a CPU version (bin/sankoff) and a GPU version (bin/cuda\_sankoff).

## Getting Started

This software have been developed and tested on Linux, and it should work on all major distributions.

### Prerequisites

You need a modern C++ compiler. In Ubuntu, you can install it by:

```
sudo apt-get install build-essential libboost-dev
```

You also need the CUDA SDK. Check how to install it on the NVIDIA page:

https://developer.nvidia.com/cuda-downloads

You also, need to Download Vienna RNA library, compile it in the folder 'ViennaRNA-2.3.3'. Check the instructions inside this folder.

### Compiling

To compile, you enter the "cuda\_sankoff" folder and type:

```
'make'
```

This works in all major Linux distributions and the 'bin/sankoff' and 'bin/cuda\_sankoff' binaries will be available.

## How to execute

Few examples:
```
#Easy test:
./bin/sankoff seqs/025.fasta

#Easy test:
./bin/cuda\_sankoff seqs/050.fasta
```

Check the 'cuda\_sankoff/seqs/' folder for other tests.

### More options

Usually the default options are the best to run. You can change and check more at:

```
./bin/sankoff -h
```

## Cite us
D. Sundfeld, G. Teodoro, J. H. Havgaard, J. Gorodkin, and A. C. M. A. Melo. Using GPU to accelerate the pairwise structural RNA alignment with base pair probabilities. *Concurrency and computation practice and experience*, Accepted, 2019. 

D. Sundfeld, J. H. Havgaard, J. Gorodkin, and A. C. M. A. Melo. CUDA-Sankoff: Using GPU to accelerate the pairwise structural RNA alignment. In *25th Euromicro International Conference on Parallel, Distributed and Network-based Processing, PDP 2017, 2017*, pages 295–302, 2017. 

## License

This project is licensed under the MIT License
