#include <iostream>
#include <omp.h>

#include "bp_probs.h"
#include "Foldalign.h"
#include "sankoff_args.h"
#include "Sankoff.h"
#include "Sequences.h"

int main(int argc, char *argv[])
{
    int error;
    int threads_num = 1;
    struct bp_prob *bp0, *bp1;

    if ((error = load_args(argc, argv, &threads_num)))
        return error;

    bp0 = new bp_prob();
    bp1 = new bp_prob();

    get_bp_prob(Sequences::get_instance()->get_seq(0), bp0);
    get_bp_prob(Sequences::get_instance()->get_seq(1), bp1);

    Sankoff fold_instance(Sequences::get_instance()->get_seq(0), Sequences::get_instance()->get_seq(1));

    omp_set_num_threads(threads_num);
    int ret = fold_instance.diag_sankoff();

    delete bp0;
    delete bp1;
    return ret;
}
