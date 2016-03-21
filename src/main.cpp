#include <iostream>
#include <omp.h>

#include "Foldalign.h"

int main(int argc, char *argv[])
{
    (void) argc;
    (void) argv;
    Foldalign fold_instance("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", "GAUCGUAGUACUGAUGAUCGUAGCAUGGUCAUGCAGUACGUUGACGUCAGUCGUUGAUGCGUACUGAGCUGUACGUCAUG", 20, 20);
    //Foldalign fold_instance("CCCCCCCCAAAAGGGGGGGG", "GUGCCAACAUUAGUUGGCAC", 20, 20);

    omp_set_num_threads(2);
    return fold_instance.diag_sankoff();
}
