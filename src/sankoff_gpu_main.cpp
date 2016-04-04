#include <iostream>
#include <omp.h>

#include "Foldalign.h"
#include "Sankoff_GPU.h"

int main(int argc, char *argv[])
{
    (void) argc;
    (void) argv;
    Sankoff_GPU fold_instance("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", "GAUCGUAGUACUGAUGAUCGUAGCAUGGUCAUGCAGUACGUUGACGUCAGUCGUUGAUGCGUACUGAGCUGUACGUCAUG");
    //Sankoff fold_instance("CCCCCCCCAAAAGGGGGGGG", "GUGCCAACAUUAGUUGGCAC");
    //Foldalign fold_instance("CCCCCCCCAAAAGGGGGGGG", "GUGCCAACAUUAGUUGGCAC", 20, 20);
    //return fold_instance.fold_align();

    //omp_set_num_threads(4);
    return fold_instance.diag_sankoff();
}
