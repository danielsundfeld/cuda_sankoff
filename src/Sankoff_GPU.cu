#include "Sankoff_GPU.h"

#include <iostream>

#include "Cost.h"

#define MAX(x, y) x > y ? x : y
Sankoff_GPU::Sankoff_GPU(const std::string &s1, const std::string &s2)
{
    int s1_l = (int) s1.length();
    int s2_l = (int) s2.length();

    //cudaMalloc(&dp_matrix, dp_matrix_calc_total_size(s1_l, s2_l) * sizeof(int));
    //cudaMalloc(&seq_ctx, sizeof(sequences));
    dp_matrix = (int*)malloc(dp_matrix_calc_total_size(s1_l, s2_l) * sizeof(int));
    seq_ctx = (sequences*)malloc(sizeof(sequences));

#define cudaMemcpy(x, y, z, w) memcpy(x, y, z)
    cudaMemcpy(&(seq_ctx->s1), s1.c_str(), (s1_l + 1) * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(&(seq_ctx->s1_l), &s1_l, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&(seq_ctx->s2), s2.c_str(), (s2_l + 1) * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(&(seq_ctx->s2_l), &s2_l, sizeof(int), cudaMemcpyHostToDevice);
}

Sankoff_GPU::~Sankoff_GPU()
{
#define cudaFree free
    cudaFree(dp_matrix);
    cudaFree(seq_ctx);
}

//! Expand one cell with position \a i, \a j, \a k, \a l
void Sankoff_GPU::expand_pos(int *dp_matrix, const int &i, const int &j, const int &k, const int &l, const sequences* const seq_ctx)
{
    int score = 0;
    const int &s1_l = seq_ctx->s1_l;
    const int &s2_l = seq_ctx->s1_l;
    const char* const &s1 = seq_ctx->s1;
    const char* const &s2 = seq_ctx->s2;

    if (j < i || j >= s1_l || l < k || l >= s2_l)
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
int Sankoff_GPU::expand_inner_matrix_diagonal_phase1(int *dp_matrix, int tid, int inner_diag, int i, int k, const sequences* const seq_ctx)
{
    int l = k + tid;
    int j = inner_diag - tid;

    expand_pos(dp_matrix, i, j, k, l, seq_ctx);
    return 0;
}

/*!
 * Expand one inner_matrix cell with coord \a i and \k, id \a tid, from one diagonal \a diag.
 * This is the first wave, from the begin to the main diagonal.
 */
int Sankoff_GPU::expand_inner_matrix_diagonal_phase2(int *dp_matrix, int tid, int inner_diag, int i, int k, const sequences* const seq_ctx)
{
    int l = inner_diag + tid;
    int j = seq_ctx->s1_l - 1 - tid;

    expand_pos(dp_matrix, i, j, k, l, seq_ctx);
    return 0;
}

/*!
 * Expand one inner_matrix cell with coord \a i and \k, id \a tid, from one diagonal \a diag.
 * This is the second wave, from the main diagonal to the end.
 */
void Sankoff_GPU::expand_inner_matrix_diag(int *dp_matrix, const int &i, const int &k, const sequences* const seq_ctx)
{
    const int &s1_l = seq_ctx->s1_l;
    const int &s2_l = seq_ctx->s2_l;

    if (i < 0)
        return;

    // First wave, from the begin to the main diagonal
    for (int inner_diag = i; inner_diag < s1_l; ++inner_diag)
    {
#pragma omp parallel for schedule(dynamic,1)
        for (int tid = 0; tid < s2_l - k; ++tid)
        {
            expand_inner_matrix_diagonal_phase1(dp_matrix, tid, inner_diag, i, k, seq_ctx);
        }
    }

    // Second wave, from the main diagonal to the end
    for (int inner_diag = k + 1; inner_diag < s2_l; ++inner_diag)
    {
#pragma omp parallel for schedule(dynamic,1)
        for (int tid = 0; tid < s1_l - inner_diag; ++tid)
        {
            expand_inner_matrix_diagonal_phase2(dp_matrix, tid, inner_diag, i, k, seq_ctx);
        }
    }
}

/*!
 * Expand one outer_matrix cell with id \a tid, from one diagonal \a diag.
 * This is the first wave, from the begin to the main diagonal.
 */
int Sankoff_GPU::expand_outer_matrix_diagonal_phase1(int *dp_matrix, int tid, int outer_diag, const sequences* const seq_ctx)
{
    int k = seq_ctx->s2_l - 1 - outer_diag + tid;
    int i = seq_ctx->s1_l - 1 - tid;

    expand_inner_matrix_diag(dp_matrix, i, k, seq_ctx);
    return 0;
}

/*!
 * Expand one outer_matrix cell with id \a tid, from one diagonal \a diag.
 * This is the second wave, from the main diagonal to the end.
 */
int Sankoff_GPU::expand_outer_matrix_diagonal_phase2(int *dp_matrix, int tid, int outer_diag, const sequences* const seq_ctx)
{
    int k = tid;
    int i = outer_diag - k;

    expand_inner_matrix_diag(dp_matrix, i, k, seq_ctx);
    return 0;
}

//! Runs Sankoff in a Diagonal way
int Sankoff_GPU::diag_sankoff()
{
    const int &s1_l = seq_ctx->s1_l;
    const int &s2_l = seq_ctx->s1_l;
    std::cout << "Sankoff_GPU:"
        << "\nseq1:\t" << seq_ctx->s1
        << "\nseq2:\t" << seq_ctx->s2
        << "\n";

    // First wave, from the begin to the main diagonal
    for (int outer_diag = 0; outer_diag <= s2_l - 1; ++outer_diag)
    {
#pragma omp parallel for schedule(dynamic,1)
        for (int tid = 0; tid <= outer_diag; ++tid)
        {
            expand_outer_matrix_diagonal_phase1(dp_matrix, tid, outer_diag, seq_ctx);
        }
    }

    // Second wave, from the main diagonal to the end
    for (int outer_diag = s1_l - 2; outer_diag >= 0 ; --outer_diag)
    {
#pragma omp parallel for schedule(dynamic,1)
        for (int tid = 0; tid <= s2_l - 1; ++tid)
        {
            expand_outer_matrix_diagonal_phase2(dp_matrix, tid, outer_diag, seq_ctx);
        }
    } //outer_diag
    //TODO copy to CPU
    std::cout << dp_matrix_get_pos(dp_matrix, 0, s1_l - 1, 0, s2_l - 1, seq_ctx) << std::endl;
    std::cout << "FIM\n";
    return 0;
}
