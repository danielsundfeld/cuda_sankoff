#include <iostream>
#include <omp.h>

#include "sankoff_args.h"
#include "Sankoff.h"
#include "Sequences.h"
#include "TimeCounter.h"

#include "dp_matrix_cell.h"

long long int dp_matrix_calc_total_size(long long int s1, long long int s2)
{
    return (((1 + s1) * s1) / 2 ) * (((1 + s2) * s2) / 2);
}

int main(int argc, char *argv[])
{
    TimeCounter t("Total execution time");
    int error;
    int threads_num = 1;
    if ((error = load_args(argc, argv, &threads_num)))
        return error;
    //Sankoff fold_instance(Sequences::get_instance()->get_seq(0), Sequences::get_instance()->get_seq(1));


	std::cout << dp_matrix_calc_total_size(3000, 3000) * sizeof(dp_matrix_cell) << std::endl;
    omp_set_num_threads(threads_num);
    //return fold_instance.diag_sankoff();
	return 0;
}
