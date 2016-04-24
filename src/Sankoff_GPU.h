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
        void expand_inner_matrix_diag(int *dp_matrix, const int &i, const int &k, const sequences* const seq_ctx);
        void expand_pos(int *dp_matrix, const int &i, const int &j, const int &k, const int &l, const sequences* const seq_ctx);

        int expand_outer_matrix_diagonal_phase1(int *dp_matrix, int tid, int outer_diag, const sequences* const seq_ctx);
        int expand_outer_matrix_diagonal_phase2(int *dp_matrix, int tid, int outer_diag, const sequences* const seq_ctx);

        int expand_inner_matrix_diagonal_phase1(int *dp_matrix, int tid, int outer_diag, int i, int k, const sequences* const seq_ctx);
        int expand_inner_matrix_diagonal_phase2(int *dp_matrix, int tid, int outer_diag, int i, int k, const sequences* const seq_ctx);

        //device_members
        int *dp_matrix;
        sequences *seq_ctx;
};
#endif
