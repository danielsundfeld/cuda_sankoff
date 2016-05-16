#include "Sankoff_GPU.h"

#include <iostream>
#include <cstring>

#include "Cost.h"

#include "DPMatrix_GPU.cu"

#define MAX(x, y) x > y ? x : y
Sankoff_GPU::Sankoff_GPU(const std::string &seq1, const std::string &seq2)
{
    h_seq_ctx = sequences();
    std::strcpy(h_seq_ctx.s1, seq1.c_str());
    std::strcpy(h_seq_ctx.s2, seq2.c_str());
    h_seq_ctx.s1_l = (int) seq1.length();
    h_seq_ctx.s2_l = (int) seq2.length();

    check_gpu_code(cudaMalloc(&dp_matrix, dp_matrix_calc_total_size(h_seq_ctx.s1_l, h_seq_ctx.s2_l) * sizeof(int)));
    check_gpu_code(cudaMalloc(&d_seq_ctx, sizeof(sequences)));
    check_gpu_code(cudaMemcpy(d_seq_ctx, &h_seq_ctx, sizeof(sequences), cudaMemcpyHostToDevice));
}

Sankoff_GPU::~Sankoff_GPU()
{
    cudaFree(dp_matrix);
    cudaFree(d_seq_ctx);
}

void Sankoff_GPU::check_gpu_code(cudaError_t code)
{
    if (code != cudaSuccess)
    {
        std::cerr << "Fatal error: " << cudaGetErrorString(code) << std::endl;
        exit(1);
    }
}

//! Expand one cell with position \a i, \a j, \a k, \a l
__device__ void sankoff_gpu_expand_pos(int *dp_matrix, const int &i, const int &j, const int &k, const int &l, sequences* seq_ctx)
{
    int score = 0;
    const char* const &s1 = seq_ctx->s1;
    const char* const &s2 = seq_ctx->s2;

    if (dp_matrix_check_border(i, j, k, l, seq_ctx) == false)
        return;

    score = MAX(score, dp_matrix_get_pos(dp_matrix, i + 1, j, k, l, seq_ctx) + Cost::gap);
    score = MAX(score, dp_matrix_get_pos(dp_matrix, i, j, k + 1, l, seq_ctx) + Cost::gap);
    score = MAX(score, dp_matrix_get_pos(dp_matrix, i, j - 1, k, l, seq_ctx) + Cost::gap);
    score = MAX(score, dp_matrix_get_pos(dp_matrix, i, j, k, l - 1, seq_ctx) + Cost::gap);
    score = MAX(score, dp_matrix_get_pos(dp_matrix, i + 1, j, k + 1, l, seq_ctx) + Cost::match_score(s1[i], s2[k]));
    score = MAX(score, dp_matrix_get_pos(dp_matrix, i, j - 1, k, l - 1, seq_ctx) + Cost::match_score(s1[j], s2[l]));
    score = MAX(score, dp_matrix_get_pos(dp_matrix, i + 1, j - 1, k, l, seq_ctx) + Cost::base_score(s1[i], s1[j]) + Cost::gap * 2);
    score = MAX(score, dp_matrix_get_pos(dp_matrix, i, j, k + 1, l - 1, seq_ctx) + Cost::base_score(s2[k], s2[l]) + Cost::gap * 2);
    score = MAX(score, dp_matrix_get_pos(dp_matrix, i + 1, j - 1, k + 1, l - 1, seq_ctx) +
            Cost::base_score(s1[i], s1[j]) + Cost::base_score(s2[k], s2[l]) +
            Cost::compensation_score(s1[i], s1[j], s2[k], s2[l]));

    for (int m = i + 1; m < j; ++m)
    {
        for (int n = k + 1; n < l; ++n)
        {
            score = MAX(score, dp_matrix_get_pos(dp_matrix, i, m, k, n, seq_ctx) + dp_matrix_get_pos(dp_matrix, m + 1, j, n + 1, l, seq_ctx));
        } //n
    } //m
    dp_matrix_put_pos(dp_matrix, i, j, k, l, score, seq_ctx);
}

//! Expand inner matrix, first wave, from the begin to the main diagonal
__device__ void sankoff_gpu_expand_inner_matrix_diagonal_phase1(int *dp_matrix, int inner_diag, int i, int k, sequences* seq_ctx)
{
    int l = k + threadIdx.x;
    int j = inner_diag - threadIdx.x;

    sankoff_gpu_expand_pos(dp_matrix, i, j, k, l, seq_ctx);
    return;
}

/*!
 * Expand one inner_matrix cell with coord \a i and \k, id \a tid, from one diagonal \a diag.
 * This is the first wave, from the begin to the main diagonal.
 */
__device__ void sankoff_gpu_expand_inner_matrix_diagonal_phase2(int *dp_matrix, int inner_diag, int i, int k, sequences* seq_ctx)
{
    int l = inner_diag + threadIdx.x;
    int j = seq_ctx->s1_l - 1 - threadIdx.x;

    sankoff_gpu_expand_pos(dp_matrix, i, j, k, l, seq_ctx);
    return;
}

/*!
 * Expand one inner_matrix cell with coord \a i and \k, id \a tid, from one diagonal \a diag.
 * This is the second wave, from the main diagonal to the end.
 */
__device__ void sankoff_gpu_expand_inner_matrix_diag(int *dp_matrix, const int &i, const int &k, sequences* seq_ctx)
{
    const int &s1_l = seq_ctx->s1_l;
    const int &s2_l = seq_ctx->s2_l;

    if (i < 0 || i >= s1_l || k < 0 || k >= s2_l)
        return;

    // First wave, from the begin to the main diagonal
    for (int inner_diag = i; inner_diag < s1_l; ++inner_diag)
    {
        sankoff_gpu_expand_inner_matrix_diagonal_phase1(dp_matrix, inner_diag, i, k, seq_ctx);
        __syncthreads();
    }

    // Second wave, from the main diagonal to the end
    for (int inner_diag = k + 1; inner_diag < s2_l; ++inner_diag)
    {
        sankoff_gpu_expand_inner_matrix_diagonal_phase2(dp_matrix, inner_diag, i, k, seq_ctx);
        __syncthreads();
    }
}

/*!
 * Expand one outer_matrix cell with id \a tid, from one diagonal \a diag.
 * This is the first wave, from the begin to the main diagonal.
 */
__global__ void sankoff_gpu_expand_outer_matrix_diagonal_phase1(int *dp_matrix, int outer_diag, sequences* seq_ctx)
{
    int k = seq_ctx->s2_l - 1 - outer_diag + blockIdx.x;
    int i = seq_ctx->s1_l - 1 - blockIdx.x;

    sankoff_gpu_expand_inner_matrix_diag(dp_matrix, i, k, seq_ctx);
    return;
}

/*!
 * Expand one outer_matrix cell with id \a tid, from one diagonal \a diag.
 * This is the second wave, from the main diagonal to the end.
 */
__global__ void sankoff_gpu_expand_outer_matrix_diagonal_phase2(int *dp_matrix, int outer_diag, sequences* seq_ctx)
{
    int k = blockIdx.x;
    int i = outer_diag - k;

    sankoff_gpu_expand_inner_matrix_diag(dp_matrix, i, k, seq_ctx);
    return;
}

//! Runs Sankoff in a Diagonal way
int Sankoff_GPU::diag_sankoff()
{
    std::cout << "Sankoff_GPU:"
        << "\nseq1: (" << h_seq_ctx.s1_l << ")\t" << h_seq_ctx.s1
        << "\nseq2: (" << h_seq_ctx.s2_l << ")\t" << h_seq_ctx.s2
        << "\n";

    int threads_num = 0;
    // First wave, from the begin to the main diagonal
    for (int outer_diag = 0; outer_diag <= h_seq_ctx.s2_l - 1; ++outer_diag)
    {
        ++threads_num;
        dim3 tn = ceil((float)threads_num/2);
        sankoff_gpu_expand_outer_matrix_diagonal_phase1<<<outer_diag + 1, tn>>>(dp_matrix, outer_diag, d_seq_ctx);
    }

    // Second wave, from the main diagonal to the end
    for (int outer_diag = h_seq_ctx.s1_l - 2; outer_diag >= 0 ; --outer_diag)
    {
        ++threads_num;
        dim3 tn = ceil((float)threads_num/2);
        sankoff_gpu_expand_outer_matrix_diagonal_phase2<<<outer_diag + 1, tn>>>(dp_matrix, outer_diag, d_seq_ctx);
    } //outer_diag
    std::cout << "Score: " << dp_matrix_get_val(dp_matrix, 0, h_seq_ctx.s1_l - 1, 0, h_seq_ctx.s2_l - 1, &h_seq_ctx) << std::endl;
    return 0;
}
