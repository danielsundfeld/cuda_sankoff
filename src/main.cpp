#include <iostream>
#include <omp.h>

#include "Foldalign.h"
#include "sankoff_args.h"
#include "Sankoff.h"
#include "Sequences.h"

int main(int argc, char *argv[])
{
    int error;
    if ((error = load_file(argc, argv)))
        return error;
    Sankoff fold_instance(Sequences::get_instance()->get_seq(0), Sequences::get_instance()->get_seq(1));

//    omp_set_num_threads(4);
    return fold_instance.diag_sankoff();
}
