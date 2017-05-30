#ifndef _SANKOFF_GPU_H
#define _SANKOFF_GPU_H
#include <string>

#include "bp_probs.h"
#include "dp_matrix_cell.h"
#include "DPMatrix_GPU.h"

class Sankoff_GPU {
    public:
        Sankoff_GPU(const std::string &seq1, const std::string &seq2);
        virtual ~Sankoff_GPU();
#ifdef __CUDACC__
        void check_gpu_code(cudaError_t code);
#endif
        int diag_sankoff(); //Run a pure sankoff algorithm

    private:
        void backtrace();
        //host members
        sequences h_seq_ctx;
        struct bp_prob *h_bp1, *h_bp2;

        //device_members
        dp_matrix_cell *dp_matrix;
        sequences *d_seq_ctx;
        struct bp_prob *d_bp1, *d_bp2;
};

#ifdef __CUDACC__
__device__ void max(dp_matrix_cell &score1, dp_matrix_cell score2, int parent);
__device__ void calculate_pos(dp_matrix_cell *dp_matrix, sequences* seq_ctx, dp_matrix_cell &score1, int i, int j, int k, int l, float extra_score, int parent);
__device__ void calculate_pos_mb(dp_matrix_cell *dp_matrix, sequences* seq_ctx, dp_matrix_cell &score1, int i, int j, int k, int l, int m, int n);
__device__ void sankoff_gpu_expand_inner_matrix_diag(dp_matrix_cell *dp_matrix, const int &i, const int &k, sequences* seq_ctx, struct bp_prob* bp1, struct bp_prob* bp2);
__device__ void sankoff_gpu_expand_pos(dp_matrix_cell *dp_matrix, const int &i, const int &j, const int &k, const int &l, sequences* seq_ctx, struct bp_prob* bp1, struct bp_prob* bp2);

__global__ void sankoff_gpu_expand_outer_matrix_diagonal_phase1(dp_matrix_cell *dp_matrix, int outer_diag, sequences* seq_ctx, struct bp_prob* bp1, struct bp_prob* bp2);
__global__ void sankoff_gpu_expand_outer_matrix_diagonal_phase2(dp_matrix_cell *dp_matrix, int outer_diag, sequences* seq_ctx, struct bp_prob* bp1, struct bp_prob* bp2);

__device__ void sankoff_gpu_expand_inner_matrix_diagonal_phase1(dp_matrix_cell *dp_matrix, int outer_diag, int i, int k, sequences* seq_ctx, struct bp_prob* bp1, struct bp_prob* bp2);
__device__ void sankoff_gpu_expand_inner_matrix_diagonal_phase2(dp_matrix_cell *dp_matrix, int outer_diag, int i, int k, sequences* seq_ctx, struct bp_prob* bp1, struct bp_prob* bp2);
#endif
#endif
