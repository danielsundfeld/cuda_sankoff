#include <iostream>
#include <omp.h>

#include "Sankoff_GPU.h"

int main(int argc, char *argv[])
{
    (void) argc;
    (void) argv;
    Sankoff_GPU fold_instance("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", "GAUCGUAGUACUGAUGAUCGUAGCAUGGUCAUGCAGUACGUUGACGUCAGUCGUUGAUGCGUACUGAGCUGUACGUCAUG");
    //Sankoff_GPU fold_instance("CCCCCCCCAAAAGGGGGGGG", "GUGCCAACAUUAGUUGGCAC");

    return fold_instance.diag_sankoff();
}
