#ifndef _SANKOFF_H
#define _SANKOFF_H
#include <string>

#include "bp_probs.h"
#include "DPMatrix.h"
#include "dp_matrix_cell.h"

class Sankoff {
    public:
        Sankoff(const std::string &seq1, const std::string &seq2);
        virtual ~Sankoff();
        int sankoff(); //Run a pure sankoff algorithm
        int diag_sankoff(); //Run a pure sankoff algorithm

        static void max(dp_matrix_cell &score1, dp_matrix_cell score2, int parent);
        void calculate_pos(dp_matrix_cell &score1, int i, int j, int k, int l, float extra_score, int parent);
        void calculate_pos_mb(dp_matrix_cell &score1, int i, int j, int k, int l, int m, int n);
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
        struct bp_prob *bp1, *bp2;

        void expand_inner_matrix(const int &i, const int &k);
        void expand_inner_matrix_diag(const int &i, const int &k);
        void expand_pos(const int &i, const int &j, const int &k, const int &l);
};
#endif
