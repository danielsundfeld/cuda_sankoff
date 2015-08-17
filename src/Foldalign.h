#ifndef _FOLDALIGN_H
#define _FOLDALIGN_H
#include <string>
#include <map>
#include <tuple>

typedef std::tuple<int, int, int, int> coord;

class Foldalign {
    public:
        Foldalign(const std::string &s1, const std::string &s2, const int lambda, const int delta);
        bool out_of_border_lambda(const int j);
        bool out_of_border_delta(const int l);
        bool out_of_border(const int i, const int j, const int k, const int l) const;
        int fold_align();
    private:
        int calculate_score(int i, int j, int k, int l);
        int calculate_mb(int i, int j, int k, int l, int m, int n);
        void print_coord(int i, int j, int k, int l) const;

        std::string m_seq1;
        int m_seq1_l;
        std::string m_seq2;
        int m_seq2_l;
        int m_lambda;
        int m_delta;

        std::map<coord, int> dp_matrix;
};
#endif
