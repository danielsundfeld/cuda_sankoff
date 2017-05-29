#ifndef _BACKTRACE_H
#define _BACKTRACE_H
#include <string>

#include "DPMatrix.h"
#include "dp_matrix_cell.h"

class DPMatrix;

class Backtrace {
    public:
        Backtrace(DPMatrix *dp_matrix, int i, int j, int k, int l, const std::string &s1, const std::string &s2);
        void run();
        void print(std::string &alignment_s1, std::string &alignment_structure, std::string &alignment_s2);

    private:
        dp_matrix_cell get_parent(const dp_matrix_cell c);
        void calculate_mb_position(float score);
        void do_backtrace_mb(int i, int j, int k, int l);

        DPMatrix *dp_matrix;
        int i, j, k, l;
        int m, n;
        const std::string s1, s2;
        std::string list_i, list_j, list_k, list_l, list_bp_left, list_bp_right;
};
#endif
