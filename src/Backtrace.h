#ifndef _BACKTRACE_H
#define _BACKTRACE_H
#include <list>
#include <string>

#include "DPMatrix.h"
#include "dp_matrix_cell.h"

class DPMatrix;

class Backtrace {
    public:
        Backtrace(DPMatrix *dp_matrix, int i, int j, int k, int l, const std::string &s1, const std::string &s2);
        void run();

    private:
        void add_last(const dp_matrix_cell c);
        dp_matrix_cell get_parent(const dp_matrix_cell c);

        DPMatrix *dp_matrix;
        int i, j, k, l;
        const std::string s1, s2;
        std::list<char> list_i, list_j, list_k, list_l, list_bp_left, list_bp_right;
};
#endif
