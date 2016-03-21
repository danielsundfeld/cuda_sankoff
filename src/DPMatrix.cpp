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
