#include <iostream>
#include <omp.h>

#include "sankoff_args.h"
#include "Sankoff.h"
#include "Sequences.h"
#include "TimeCounter.h"

int main(int argc, char *argv[])
{
    TimeCounter t("Total execution time");
    int error;
    int threads_num = 1;
    if ((error = load_args(argc, argv, &threads_num)))
        return error;
    Sankoff fold_instance(Sequences::get_instance()->get_seq(0), Sequences::get_instance()->get_seq(1));

    omp_set_num_threads(threads_num);
    return fold_instance.diag_sankoff();
}
