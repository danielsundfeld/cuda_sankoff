#ifndef _SANKOFF_GPU_H
#define _SANKOFF_GPU_H
#include <string>

#include "DPMatrix.h"

class Sankoff_GPU {
    public:
        Sankoff_GPU(const std::string &seq1, const std::string &seq2);
        virtual ~Sankoff_GPU();
        int sankoff(); //Run a pure sankoff algorithm
        int diag_sankoff(); //Run a pure sankoff algorithm

        void print_score_dep(int i, int j, int k, int l) const;
        void print_mb_dep(int i, int j, int k, int l, int m, int n) const;
        void print_orig(int i, int j, int k, int l) const;
        void print_index(int i, int j, int k, int l) const;

    protected:
        DPMatrix dp_matrix;
        std::string s1;
        int s1_l;
        std::string s2;
        int s2_l;

    private:
        void expand_inner_matrix(const int &i, const int &k);
        void expand_inner_matrix_diag(const int &i, const int &k);
        void expand_pos(const int &i, const int &j, const int &k, const int &l);
};
#endif
