#ifndef _FOLDALIGN_H
#define _FOLDALIGN_H
#include <string>
#include <map>
#include <tuple>

typedef std::tuple<int, int, int, int> index;

class Foldalign {
    public:
        Foldalign(const std::string &seq1, const std::string &seq2, const int lmbd, const int dlt);
        int sankoff(); //Run a pure sankoff algorithm
        int diag_sankoff(); //Run a pure sankoff algorithm
        int fold_align(); //Run with the Foldalign heuristic

    private:
        void expand_inner_matrix(const int &i, const int &k);
        void expand_inner_matrix_diag(const int &i, const int &k);
        void expand_pos(const int &i, const int &j, const int &k, const int &l);
        void print_score_dep(int i, int j, int k, int l) const;
        void print_mb_dep(int i, int j, int k, int l, int m, int n) const;
        void print_orig(int i, int j, int k, int l) const;
        void print_index(int i, int j, int k, int l) const;

        std::string s1;
        int s1_l;
        std::string s2;
        int s2_l;
        int lambda;
        int delta;

        std::map<index, int> dp_matrix;
};
#endif
