#include "sankoff_args.h"

#include <iostream>
#include <stdlib.h>
#include <unistd.h>

#include "read_fasta.h"

int usage(char name[])
{
    std::cerr << "Usage: " << name << "[-t threads_num] <arq.fasta>\n";
    return -1;
}

//! Load sequences to the Sequences singleton, return error if necessary
int load_args(int argc, char *argv[], int *threads_num)
{
    int opt;
    while ((opt = getopt(argc, argv, "t:")) != -1)
    {
        switch (opt)
        {
            case 't':
                if (threads_num)
                    *threads_num = atoi(optarg);
                break;
            default:
                return usage(argv[0]);
        }
    }
    if (optind >= argc)
        return usage(argv[0]);
    return read_fasta_file(argv[optind], 2);
}
