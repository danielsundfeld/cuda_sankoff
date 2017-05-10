#ifndef _DPMATRIX_H
#define _DPMATRIX_H
#include <string>

#include "Backtrace.h"
#include "dp_matrix_cell.h"

class DPMatrix {
    public:
        DPMatrix(const int &s1_l, const int &s2_l);
        ~DPMatrix();
        dp_matrix_cell get_pos(const int &i, const int &j, const int &k, const int &l) const;
        void put_pos(const int &i, const int &j, const int &k, const int &l, const dp_matrix_cell &val);
        long long int get_total_size() const;
        bool check_border(const int &i, const int &j, const int &k, const int &l) const;
        void backtrace(const std::string &s1, const std::string &s2);
        dp_matrix_cell* get_dp_matrix() { return dp_matrix; };

    private:
        int calc_delta(int i, int j, int k, int l) const;
        long long int calc_total_size(long long int s1, long long int s2) const;
        dp_matrix_cell get_parent(const dp_matrix_cell c, int &i, int &j, int &k, int &l, const std::string &s1, const std::string &s2);
        void add_last(const dp_matrix_cell c, int &i, int &k, const std::string &s1, const std::string &s2);

        dp_matrix_cell *dp_matrix;
        const int s1_l;
        const int s2_l;
};
#endif
