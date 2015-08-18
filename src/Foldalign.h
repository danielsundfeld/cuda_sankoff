#ifndef _FOLDALIGN_H
#define _FOLDALIGN_H
#include <string>
#include <map>
#include <tuple>

typedef std::tuple<int, int, int, int> index;

class Foldalign {
    public:
        Foldalign(const std::string &s1, const std::string &s2, const int lambda, const int delta);
        bool out_of_border_lambda(const int j) const;
        bool out_of_border_delta(const int l) const;
        bool out_of_border(const int i, const int j, const int k, const int l) const;
        int fold_align();
    private:
        void print_score_dep(int i, int j, int k, int l) const;
        void print_mb_dep(int i, int j, int k, int l, int m, int n) const;
        void print_orig(int i, int j, int k, int l) const;
        void print_index(int i, int j, int k, int l) const;

        std::string m_seq1;
        int m_seq1_l;
        std::string m_seq2;
        int m_seq2_l;
        int m_lambda;
        int m_delta;

        std::map<index, int> dp_matrix;
};
#endif
