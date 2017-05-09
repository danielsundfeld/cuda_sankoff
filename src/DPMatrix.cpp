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

dp_matrix_cell DPMatrix::get_parent(const dp_matrix_cell c, int &i, int &j, int &k, int &l, const std::string &s1, const std::string &s2)
{
    switch (c.parent)
    {
        case GapI:
            list_i.push_back(s1[i]);
            list_bp_left.push_back('.');
            list_k.push_back('-');
            i += 1;
            break;

        case GapK:
            list_k.push_back(s2[k]);
            list_bp_left.push_back('.');
            list_i.push_back('-');
            k += 1;
            break;

        case GapJ:
            list_j.push_front(s1[j]);
            list_bp_right.push_front('.');
            list_l.push_front('-');
            j -= 1;
            break;

        case GapL:
            list_l.push_front(s2[l]);
            list_bp_right.push_front('.');
            list_j.push_front('-');
            l -= 1;
            break;

        case UnpairedIK:
            list_i.push_back(s1[i]);
            list_bp_left.push_back('.');
            list_k.push_back(s2[k]);
            i += 1;
            k += 1;
            break;

        case UnpairedJL:
            list_j.push_front(s1[j]);
            list_bp_right.push_front('.');
            list_l.push_front(s2[l]);
            j -= 1;
            l -= 1;
            break;

        case PairedGapS1:
            list_i.push_back(s1[i]);
            list_j.push_front(s1[j]);
            list_bp_left.push_back('(');
            list_bp_right.push_front(')');
            list_k.push_back('-');
            list_l.push_front('-');
            i += 1;
            j -= 1;
            break;

        case PairedGapS2:
            list_k.push_back(s2[k]);
            list_l.push_front(s2[l]);
            list_bp_left.push_back('(');
            list_bp_right.push_front(')');
            list_i.push_back('-');
            list_j.push_front('-');
            k += 1;
            l -= 1;
            break;

        case Paired:
            list_i.push_back(s1[i]);
            list_j.push_front(s1[j]);
            list_bp_left.push_back('(');
            list_bp_right.push_front(')');
            list_k.push_back(s2[k]);
            list_l.push_front(s2[l]);
            i += 1;
            j -= 1;
            k += 1;
            l -= 1;
            break;
    }
    return get_pos(i, j, k, l);
}

void DPMatrix::add_last(const dp_matrix_cell c, int &i, int &k, const std::string &s1, const std::string &s2)
{
    if (c.parent != NullParent)
        return;

    list_i.push_back(s1[i]);
    list_k.push_back(s2[k]);
    list_bp_left.push_back('.');
}

void print_list(const std::list<char> &list)
{
    for (std::list<char>::const_iterator it = list.begin(); it != list.end(); ++it)
        std::cout << *it;
}

void DPMatrix::backtrace(const std::string &s1, const std::string &s2)
{
    int i = 0;
    int j = s1_l - 1;
    int k = 0;
    int l = s2_l - 1;

    dp_matrix_cell c = get_pos(i, j, k, l);

    while (check_border(i, j, k, l) && c.parent != NullParent)
    {
        printf("%f (%s) - %d %d %d %d\n", c.score, parent_str[(int)c.parent], i, k, j, l);
        c = get_parent(c, i, j, k, l, s1, s2);
    }
    add_last(c, i, k, s1, s2);
    printf("Fim: %f (%s) - %d %d %d %d\n", c.score, parent_str[(int)c.parent], i, k, j, l);

    print_list(list_i);
    print_list(list_j);
    std::cout << std::endl;
    print_list(list_bp_left);
    print_list(list_bp_right);
    std::cout << std::endl;
    print_list(list_k);
    print_list(list_l);
    std::cout << std::endl;
}
