#include <iostream>

#include "sankoff_args.h"
#include "Sankoff_GPU.h"
#include "Sequences.h"
#include "TimeCounter.h"

int main(int argc, char *argv[])
{
    TimeCounter t("Total execution time");
    int error;
    if ((error = load_args(argc, argv)))
        return error;
    Sankoff_GPU fold_instance(Sequences::get_instance()->get_seq(0), Sequences::get_instance()->get_seq(1));

    return fold_instance.diag_sankoff();
}
