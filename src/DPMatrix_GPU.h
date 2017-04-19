#ifndef _DPMATRIX_GPU_H
#define _DPMATRIX_GPU_H
#include "dp_matrix_cell.h"

#define MAX_SEQ_SIZE 1024

struct sequences {
    char s1[MAX_SEQ_SIZE];
    int s1_l;
    char s2[MAX_SEQ_SIZE];
    int s2_l;
};

long long int dp_matrix_calc_total_size(long long int s1, long long int s2);
#ifdef __CUDACC__
__device__ __host__ int dp_matrix_calc_delta(int i, int j, int k, int l, sequences* seq_ctx);
__device__ __host__ dp_matrix_cell dp_matrix_get_pos(dp_matrix_cell *dp_matrix, const int &i, const int &j, const int &k, const int &l, sequences* seq_ctx);
__device__ __host__ void dp_matrix_put_pos(dp_matrix_cell *dp_matrix, const int &i, const int &j, const int &k, const int &l, const float &val, sequences* seq_ctx);
__host__ float dp_matrix_get_val(dp_matrix_cell *dp_matrix, const int &i, const int &j, const int &k, const int &l, sequences* seq_ctx);

__device__ __host__ inline bool dp_matrix_check_border(const int &i, const int &j, const int &k, const int &l, sequences* seq_ctx)
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
#endif
#endif
