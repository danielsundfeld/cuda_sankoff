#include <iostream>
#include <omp.h>

#include "Foldalign.h"
#include "Sankoff.h"

int main(int argc, char *argv[])
{
    (void) argc;
    (void) argv;
    Sankoff fold_instance("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", "GAUCGUAGUACUGAUGAUCGUAGCAUGGUCAUGCAGUACGUUGACGUCAGUCGUUGAUGCGUACUGAGCUGUACGUCAUG");
    //Sankoff fold_instance("CCCCCCCCAAAAGGGGGGGG", "GUGCCAACAUUAGUUGGCAC");
    //Foldalign fold_instance("CCCCCCCCAAAAGGGGGGGG", "GUGCCAACAUUAGUUGGCAC", 20, 20);
    //return fold_instance.fold_align();

    omp_set_num_threads(4);
    return fold_instance.diag_sankoff();
}
