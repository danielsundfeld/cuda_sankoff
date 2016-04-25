#ifndef _SANKOFF_GPU_H
#define _SANKOFF_GPU_H
#include <string>

#include "DPMatrix_GPU.h"

class Sankoff_GPU {
    public:
        Sankoff_GPU(const std::string &seq1, const std::string &seq2);
        virtual ~Sankoff_GPU();
        int diag_sankoff(); //Run a pure sankoff algorithm

    private:
        //host members
        std::string s1;
        std::string s2;
        int s1_l;
        int s2_l;

        //device_members
        int *dp_matrix;
        sequences *seq_ctx;
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
