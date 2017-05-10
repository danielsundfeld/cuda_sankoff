#include "DPMatrix.h"

#include <iostream>
#include <stdio.h>

#include "dp_matrix_cell.h"

DPMatrix::DPMatrix(const int &s1_l, const int &s2_l)
: s1_l(s1_l),
  s2_l(s2_l)
{
    dp_matrix = new dp_matrix_cell[calc_total_size(s1_l, s2_l)]();
}

DPMatrix::~DPMatrix()
{
    delete[] dp_matrix;
}

long long int DPMatrix::get_total_size() const
{
    return calc_total_size(s1_l, s2_l);
}

dp_matrix_cell DPMatrix::get_pos(const int &i, const int &j, const int &k, const int &l) const
{
    if (check_border(i, j, k, l) == false)
    {
        dp_matrix_cell ret = dp_matrix_cell();
        ret.score = -1024;
        return ret;
    }

    return dp_matrix[calc_delta(i, j, k, l)];
}

void DPMatrix::put_pos(const int &i, const int &j, const int &k, const int &l, const dp_matrix_cell &val)
{
    if (check_border(i, j, k, l) == false)
        return;

    dp_matrix[calc_delta(i, j, k, l)] = val;
}

int DPMatrix::calc_delta(int i, int j, int k, int l) const
{
    i = s1_l - i;
    j = s1_l - j;
    k = s2_l - k;
    l = s2_l - l;

    int delta_i = ((1 + s2_l) * s2_l / 2) * (i * (i - 1) / 2);
    int delta_k = i * (k * (k - 1) / 2);
    int delta_mi = (j - 1) * k + l - 1;
    //std::cout << delta_i + delta_k + delta_mi << " (" << delta_i << " " << delta_k << " " << delta_mi << ") " << i << " " << j << " " << k << " " << l << "\n";
    return delta_i + delta_k + delta_mi;
}

long long int DPMatrix::calc_total_size(long long int s1, long long int s2) const
{
    return (((1 + s1) * s1) / 2 ) * (((1 + s2) * s2) / 2);
}

bool DPMatrix::check_border(const int &i, const int &j, const int &k, const int &l) const
{
    if (j < i)
        return false;
    if (l < k)
        return false;
    if (i < 0 || j < 0 || k < 0 || l < 0)
        return false;
    if (i >= s1_l || j >= s1_l || k >= s2_l || l >= s2_l)
        return false;
    return true;
}

void DPMatrix::backtrace(const std::string &s1, const std::string &s2)
{
    std::string alignment_s1, alignment_structure, alignment_s2;
    Backtrace bc(this, 0, s1_l - 1, 0, s2_l - 1, s1, s2);
    bc.run();
    bc.print(alignment_s1, alignment_structure, alignment_s2);
    std::cout << alignment_s1 << "\n" << alignment_structure << "\n" << alignment_s2 << "\n";
}
