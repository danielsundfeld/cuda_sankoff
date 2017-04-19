#include "DPMatrix_GPU.h"

#include "dp_matrix_cell.h"

long long int dp_matrix_calc_total_size(long long int s1, long long int s2)
{
    return (((1 + s1) * s1) / 2 ) * (((1 + s2) * s2) / 2);
}

__device__ __host__ int dp_matrix_calc_delta(int i, int j, int k, int l, sequences* seq_ctx)
{
    const int &s1_l = seq_ctx->s1_l;
    const int &s2_l = seq_ctx->s2_l;

    i = s1_l - i;
    j = s1_l - j;
    k = s2_l - k;
    l = s2_l - l;

    int delta_i = ((1 + s2_l) * s2_l / 2) * (i * (i - 1) / 2);
    int delta_k = i * (k * (k - 1) / 2);
    int delta_mi = (j - 1) * k + l - 1;
    return delta_i + delta_k + delta_mi;
}

__device__ __host__ dp_matrix_cell dp_matrix_get_pos(dp_matrix_cell *dp_matrix, const int &i, const int &j, const int &k, const int &l, sequences* seq_ctx)
{
    if (dp_matrix_check_border(i, j, k, l, seq_ctx) == false)
    {
        dp_matrix_cell ret = dp_matrix_cell();
        ret.score = -1024;
        return ret;
    }

    return dp_matrix[dp_matrix_calc_delta(i, j, k, l, seq_ctx)];
}

__device__ __host__ void dp_matrix_put_pos(dp_matrix_cell *dp_matrix, const int &i, const int &j, const int &k, const int &l, const dp_matrix_cell &val, sequences* seq_ctx)
{
    dp_matrix[dp_matrix_calc_delta(i, j, k, l, seq_ctx)] = val;
}

__host__ float dp_matrix_get_val(dp_matrix_cell *dp_matrix, const int &i, const int &j, const int &k, const int &l, sequences* seq_ctx)
{
    if (dp_matrix_check_border(i, j, k, l, seq_ctx) == false)
        return -1024;
    dp_matrix_cell val = {};
    cudaMemcpy((void*)(&val), (void*)(&dp_matrix[dp_matrix_calc_delta(i, j, k, l, seq_ctx)]), sizeof(dp_matrix_cell), cudaMemcpyDeviceToHost);
    return val.score;
}
