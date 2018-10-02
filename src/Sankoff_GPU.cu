#include "Sankoff_GPU.h"

#include <iostream>
#include <cstring>

#include "bp_probs.h"
#include "Cost.h"
#include "dp_matrix_cell.h"
#include "DPMatrix.h"
#include "DPMatrix_GPU.cu"
#include "TimeCounter.h"

#define MAX(x, y) x > y ? x : y
Sankoff_GPU::Sankoff_GPU(const std::string &seq1, const std::string &seq2)
{
    h_seq_ctx = sequences();
    std::strcpy(h_seq_ctx.s1, seq1.c_str());
    std::strcpy(h_seq_ctx.s2, seq2.c_str());
    h_seq_ctx.s1_l = (int) seq1.length();
    h_seq_ctx.s2_l = (int) seq2.length();

    h_bp1 = new bp_prob();
    h_bp2 = new bp_prob();

    get_bp_prob(seq1, h_bp1);
    get_bp_prob(seq2, h_bp2);

    size_t dp_matrix_size = dp_matrix_calc_total_size(h_seq_ctx.s1_l, h_seq_ctx.s2_l) * sizeof(dp_matrix_cell);
    check_gpu_code(cudaMalloc(&dp_matrix, dp_matrix_size));
    check_gpu_code(cudaMemset(dp_matrix, 0, dp_matrix_size));
    check_gpu_code(cudaMalloc(&d_seq_ctx, sizeof(sequences)));
    check_gpu_code(cudaMemcpy(d_seq_ctx, &h_seq_ctx, sizeof(sequences), cudaMemcpyHostToDevice));

    check_gpu_code(cudaMalloc(&d_bp1, sizeof(struct bp_prob)));
    check_gpu_code(cudaMemcpy(d_bp1, h_bp1, sizeof(struct bp_prob), cudaMemcpyHostToDevice));
    check_gpu_code(cudaMalloc(&d_bp2, sizeof(struct bp_prob)));
    check_gpu_code(cudaMemcpy(d_bp2, h_bp2, sizeof(struct bp_prob), cudaMemcpyHostToDevice));
}

Sankoff_GPU::~Sankoff_GPU()
{
    delete h_bp1;
    delete h_bp2;
    cudaFree(dp_matrix);
    cudaFree(d_seq_ctx);
    cudaFree(d_bp1);
    cudaFree(d_bp2);
}

void Sankoff_GPU::check_gpu_code(cudaError_t code)
{
    if (code != cudaSuccess)
    {
        std::cerr << "Fatal error: " << cudaGetErrorString(code) << std::endl;
        exit(1);
    }
}

void copy_dp_matrix_from_gpu(dp_matrix_cell *h_dp_matrix, dp_matrix_cell *d_dp_matrix, long long int size)
{
    cudaError_t code = cudaMemcpy(h_dp_matrix, d_dp_matrix, sizeof(dp_matrix_cell) * size, cudaMemcpyDeviceToHost);
    if (code != cudaSuccess)
    {
        std::cerr << "Fatal error copying the dp_matrix from gpu: " << cudaGetErrorString(code) << std::endl;
        exit(1);
    }
}

__device__ void max(dp_matrix_cell &score1, dp_matrix_cell score2, int parent)
{
    if (score2.score > score1.score || (score1.parent == NullParent))
    {
        score1.score = score2.score;
        score1.parent = parent;
    }
}

__device__ void calculate_pos(dp_matrix_cell *dp_matrix, sequences* seq_ctx, dp_matrix_cell &score1, int i, int j, int k, int l, float extra_score, int parent)
{
    if (dp_matrix_check_border(i, j, k, l, seq_ctx) == false)
        return;

    dp_matrix_cell score2(dp_matrix_get_pos(dp_matrix, i, j, k, l, seq_ctx));
    score2.score += extra_score;
    max(score1, score2, parent);
}

__device__ void calculate_pos_mb(dp_matrix_cell *dp_matrix, sequences* seq_ctx, dp_matrix_cell &score1, int i, int j, int k, int l, int m, int n)
{
    dp_matrix_cell mb_right(dp_matrix_get_pos(dp_matrix, m + 1, j, n + 1, l, seq_ctx));
    if (mb_right.parent != Paired)
        return;

    dp_matrix_cell mb_left(dp_matrix_get_pos(dp_matrix, i, m, k, n, seq_ctx));
    dp_matrix_cell mb;
    mb.score = mb_left.score + mb_right.score;
    mb.parent = Multibranch;

    max(score1, mb, Multibranch);
}

//! Expand one cell with position \a i, \a j, \a k, \a l
__device__ void sankoff_gpu_expand_pos(dp_matrix_cell *dp_matrix, const int &i, const int &j, const int &k, const int &l, sequences* seq_ctx, struct bp_prob* bp1, struct bp_prob* bp2)
{
    dp_matrix_cell score = dp_matrix_cell();
    const char* const &s1 = seq_ctx->s1;
    const char* const &s2 = seq_ctx->s2;

    if (dp_matrix_check_border(i, j, k, l, seq_ctx) == false)
        return;

    float s1_score = bp1->m[i+1][j+1];
    float s2_score = bp2->m[k+1][l+1];
    if (s1_score == 0 && s2_score == 0)
    {
        s1_score = -1024;
        s2_score = -1024;
    }

    calculate_pos(dp_matrix, seq_ctx, score, i + 1, j, k, l, Cost::gap, GapI);
    calculate_pos(dp_matrix, seq_ctx, score, i, j, k + 1, l, Cost::gap, GapK);
    calculate_pos(dp_matrix, seq_ctx, score, i, j - 1, k, l, Cost::gap, GapJ);
    calculate_pos(dp_matrix, seq_ctx, score, i, j, k, l - 1, Cost::gap, GapL);
    calculate_pos(dp_matrix, seq_ctx, score, i + 1, j, k + 1, l, Cost::match_score(s1[i], s2[k]), UnpairedIK);
    calculate_pos(dp_matrix, seq_ctx, score, i, j - 1, k, l - 1, Cost::match_score(s1[j], s2[l]), UnpairedJL);
    calculate_pos(dp_matrix, seq_ctx, score, i + 1, j - 1, k, l, s1_score + Cost::gap * 2, PairedGapS1);
    calculate_pos(dp_matrix, seq_ctx, score, i, j, k + 1, l - 1, s2_score + Cost::gap * 2, PairedGapS2);
    calculate_pos(dp_matrix, seq_ctx, score, i + 1, j - 1, k + 1, l - 1, s1_score + s2_score +
            Cost::compensation_score(s1[i], s1[j], s2[k], s2[l]), Paired);

    for (int m = i + 1; m < j; ++m)
    {
        for (int n = k + 1; n < l; ++n)
        {
            calculate_pos_mb(dp_matrix, seq_ctx, score, i, j, k, l, m, n);
        } //n
    } //m

    dp_matrix_put_pos(dp_matrix, i, j, k, l, score, seq_ctx);
}

//! Expand inner matrix, first wave, from the begin to the main diagonal
__device__ void sankoff_gpu_expand_inner_matrix_diagonal_phase1(dp_matrix_cell *dp_matrix, int inner_diag, int i, int k, sequences* seq_ctx, struct bp_prob* bp1, struct bp_prob* bp2)
{
    int l = k + threadIdx.x;
    int j = inner_diag - threadIdx.x;

    sankoff_gpu_expand_pos(dp_matrix, i, j, k, l, seq_ctx, bp1, bp2);
    return;
}

/*!
 * Expand one inner_matrix cell with coord \a i and \k, id \a tid, from one diagonal \a diag.
 * This is the first wave, from the begin to the main diagonal.
 */
__device__ void sankoff_gpu_expand_inner_matrix_diagonal_phase2(dp_matrix_cell *dp_matrix, int inner_diag, int i, int k, sequences* seq_ctx, struct bp_prob* bp1, struct bp_prob* bp2)
{
    int l = inner_diag + threadIdx.x;
    int j = seq_ctx->s1_l - 1 - threadIdx.x;

    sankoff_gpu_expand_pos(dp_matrix, i, j, k, l, seq_ctx, bp1, bp2);
    return;
}

/*!
 * Expand one inner_matrix cell with coord \a i and \k, id \a tid, from one diagonal \a diag.
 * This is the second wave, from the main diagonal to the end.
 */
__device__ void sankoff_gpu_expand_inner_matrix_diag(dp_matrix_cell *dp_matrix, const int &i, const int &k, sequences* seq_ctx, struct bp_prob* bp1, struct bp_prob* bp2)
{
    const int &s1_l = seq_ctx->s1_l;
    const int &s2_l = seq_ctx->s2_l;

    if (i < 0 || i >= s1_l || k < 0 || k >= s2_l)
        return;

    // First wave, from the begin to the main diagonal
    for (int inner_diag = i; inner_diag < s1_l; ++inner_diag)
    {
        sankoff_gpu_expand_inner_matrix_diagonal_phase1(dp_matrix, inner_diag, i, k, seq_ctx, bp1, bp2);
        __syncthreads();
    }

    // Second wave, from the main diagonal to the end
    for (int inner_diag = k + 1; inner_diag < s2_l; ++inner_diag)
    {
        sankoff_gpu_expand_inner_matrix_diagonal_phase2(dp_matrix, inner_diag, i, k, seq_ctx, bp1, bp2);
        __syncthreads();
    }
}

/*!
 * Expand one outer_matrix cell with id \a tid, from one diagonal \a diag.
 * This is the first wave, from the begin to the main diagonal.
 */
__global__ void sankoff_gpu_expand_outer_matrix_diagonal_phase1(dp_matrix_cell *dp_matrix, int outer_diag, sequences* seq_ctx, struct bp_prob* bp1, struct bp_prob* bp2)
{
    int k = seq_ctx->s2_l - 1 - outer_diag + blockIdx.x;
    int i = seq_ctx->s1_l - 1 - blockIdx.x;

    sankoff_gpu_expand_inner_matrix_diag(dp_matrix, i, k, seq_ctx, bp1, bp2);
    return;
}

/*!
 * Expand one outer_matrix cell with id \a tid, from one diagonal \a diag.
 * This is the second wave, from the main diagonal to the end.
 */
__global__ void sankoff_gpu_expand_outer_matrix_diagonal_phase2(dp_matrix_cell *dp_matrix, int outer_diag, sequences* seq_ctx, struct bp_prob* bp1, struct bp_prob* bp2)
{
    int k = blockIdx.x;
    int i = outer_diag - k;

    sankoff_gpu_expand_inner_matrix_diag(dp_matrix, i, k, seq_ctx, bp1, bp2);
    return;
}

void Sankoff_GPU::backtrace()
{
    TimeCounter tb("Backtrace total time");
    TimeCounter *t2 = new TimeCounter("Initializing CPU memory");
    DPMatrix h_dp_matrix(h_seq_ctx.s1_l, h_seq_ctx.s2_l);
    delete t2;
    TimeCounter *t1 = new TimeCounter("Copy from GPU");
    copy_dp_matrix_from_gpu(h_dp_matrix.get_dp_matrix(), dp_matrix, h_dp_matrix.get_total_size());
    delete t1;
    h_dp_matrix.backtrace(h_seq_ctx.s1, h_seq_ctx.s2);
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
        sankoff_gpu_expand_outer_matrix_diagonal_phase1<<<outer_diag + 1, tn>>>(dp_matrix, outer_diag, d_seq_ctx, d_bp1, d_bp2);
        //cudaDeviceSynchronize();
    }

    // Second wave, from the main diagonal to the end
    for (int outer_diag = h_seq_ctx.s1_l - 2; outer_diag >= 0 ; --outer_diag)
    {
        ++threads_num;
        dim3 tn = ceil((float)threads_num/2);
        sankoff_gpu_expand_outer_matrix_diagonal_phase2<<<outer_diag + 1, tn>>>(dp_matrix, outer_diag, d_seq_ctx, d_bp1, d_bp2);
        //cudaDeviceSynchronize();
    } //outer_diag
    cudaDeviceSynchronize();
    std::cout << "Score: " << dp_matrix_get_val(dp_matrix, 0, h_seq_ctx.s1_l - 1, 0, h_seq_ctx.s2_l - 1, &h_seq_ctx) << std::endl;

    backtrace();
    return 0;
}
