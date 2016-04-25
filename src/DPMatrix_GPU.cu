#include "DPMatrix_GPU.h"

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

__device__ __host__ bool dp_matrix_check_border(const int &i, const int &j, const int &k, const int &l, sequences* seq_ctx)
{
    const int &s1_l = seq_ctx->s1_l;
    const int &s2_l = seq_ctx->s2_l;

    if (j < i)
        return false;
    if (l < k)
        return false;
    if (i < 0 || j < 0 || k < 0 || l < 0)
        return false;
    if (i >= s1_l || j >= s1_l || k >= s2_l || l >= s2_l)
        return false;
    return true;
}

__device__ __host__ int dp_matrix_get_pos(int *dp_matrix, const int &i, const int &j, const int &k, const int &l, sequences* seq_ctx)
{
    if (dp_matrix_check_border(i, j, k, l, seq_ctx) == false)
        return -1024;

    return dp_matrix[dp_matrix_calc_delta(i, j, k, l, seq_ctx)];
}

__device__ __host__ void dp_matrix_put_pos(int *dp_matrix, const int &i, const int &j, const int &k, const int &l, const int &val, sequences* seq_ctx)
{
    dp_matrix[dp_matrix_calc_delta(i, j, k, l, seq_ctx)] = val;
}

__host__ int dp_matrix_get_val(int *dp_matrix, const int &i, const int &j, const int &k, const int &l, sequences* sq)
{
    sequences seq_ctx;
    //TODO: optimize
    cudaMemcpy((void*)&seq_ctx, sq, sizeof(sequences), cudaMemcpyDeviceToHost);
    if (dp_matrix_check_border(i, j, k, l, &seq_ctx) == false)
        return -1024;
    int val = 0;
    cudaMemcpy((void*)(&val), (void*)(&dp_matrix[dp_matrix_calc_delta(i, j, k, l, &seq_ctx)]), sizeof(int), cudaMemcpyDeviceToHost);
    return val;
}
