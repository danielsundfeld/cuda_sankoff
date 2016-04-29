#include "sankoff_args.h"

#include <iostream>

#include "read_fasta.h"

int usage(char name[])
{
    std::cerr << "Usage: " << name << " <arq.fasta>\n";
    return -1;
}

//! Load sequences to the Sequences singleton, return error if necessary
int load_file(int argc, char *argv[])
{
    if (argc == 1)
        return usage(argv[0]);
    return read_fasta_file(argv[1]);
}
