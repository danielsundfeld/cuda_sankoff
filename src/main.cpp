#include <iostream>
#include <omp.h>

#include "Sankoff.h"

int main(int argc, char *argv[])
{
    (void) argc;
    (void) argv;
    Sankoff fold_instance("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", "GAUCGUAGUACUGAUGAUCGUAGCAUGGUCAUGCAGUACGUUGACGUCAGUCGUUGAUGCGUACUGAGCUGUACGUCAUG", 20, 20);
    //Sankoff fold_instance("CCCCCCCCAAAAGGGGGGGG", "GUGCCAACAUUAGUUGGCAC", 20, 20);

    //omp_set_num_threads(4);
    return fold_instance.diag_sankoff();
}
