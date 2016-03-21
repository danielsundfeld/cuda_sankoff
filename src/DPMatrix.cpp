#include "DPMatrix.h"

DPMatrix::DPMatrix(const int &s1_l, const int &s2_l)
: s1_l(s1_l),
  s2_l(s2_l)
{

}

int DPMatrix::get_pos(const int &i, const int &j, const int &k, const int &l)
{
    return matrix[index(i, j, k, l)];
}

void DPMatrix::put_pos(const int &i, const int &j, const int &k, const int &l, const int &val)
{
    matrix[index(i, j, k, l)] = val;
}

int DPMatrix::calc_delta(int i, int j, int k, int l) const
{
    i = s1_l - i;
    j = s1_l - j;
    k = s2_l - k;
    l = s2_l - l;

    int delta_i = ((1 + s1_l) * s1_l / 2) * (i * (i - 1) / 2);
    int delta_k = i * (k * (k - 1) / 2);
    int delta_mi = (j - 1) * k + l - 1;
    //std::cout << delta_i + delta_k + delta_mi << " (" << delta_i << " " << delta_k << " " << delta_mi << ") " << i << " " << j << " " << k << " " << l << "\n";
    return delta_i + delta_k + delta_mi;
}
