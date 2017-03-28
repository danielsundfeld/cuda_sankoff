#ifndef _SANKOFF_GPU_H
#define _SANKOFF_GPU_H
#include <string>

#include "bp_probs.h"
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
        //host members
        sequences h_seq_ctx;
        struct bp_prob *h_bp1, *h_bp2;

        //device_members
        float *dp_matrix;
        sequences *d_seq_ctx;
        struct bp_prob *d_bp1, *d_bp2;
};

#ifdef __CUDACC__
__device__ void sankoff_gpu_expand_inner_matrix_diag(int *dp_matrix, const int &i, const int &k, sequences* seq_ctx);
__device__ void sankoff_gpu_expand_pos(int *dp_matrix, const int &i, const int &j, const int &k, const int &l, sequences* seq_ctx);

__global__ void sankoff_gpu_expand_outer_matrix_diagonal_phase1(int *dp_matrix, int outer_diag, sequences* seq_ctx);
__global__ void sankoff_gpu_expand_outer_matrix_diagonal_phase2(int *dp_matrix, int outer_diag, sequences* seq_ctx);

__device__ void sankoff_gpu_expand_inner_matrix_diagonal_phase1(int *dp_matrix, int outer_diag, int i, int k, sequences* seq_ctx);
__device__ void sankoff_gpu_expand_inner_matrix_diagonal_phase2(int *dp_matrix, int outer_diag, int i, int k, sequences* seq_ctx);
#endif
#endif
