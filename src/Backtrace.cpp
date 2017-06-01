#include "Backtrace.h"

#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>

Backtrace::Backtrace(DPMatrix *dp_matrix, int i, int j, int k, int l, const std::string &s1, const std::string &s2)
: dp_matrix(dp_matrix),
  i(i),
  j(j),
  k(k),
  l(l),
  s1(s1),
  s2(s2)
{

}

dp_matrix_cell Backtrace::get_parent(const dp_matrix_cell c)
{
    //small reminder: insert(0, 1, c) is push_front: insert at position 0, one char 'c'
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
            list_j.insert(0, 1, s1[j]);
            list_bp_right.insert(0, 1, '.');
            list_l.insert(0, 1, '-');
            j -= 1;
            break;

        case GapL:
            list_l.insert(0, 1, s2[l]);
            list_bp_right.insert(0, 1, '.');
            list_j.insert(0, 1, '-');
            l -= 1;
            break;

        case NullParent:
        case UnpairedIK:
            list_i.push_back(s1[i]);
            list_bp_left.push_back('.');
            list_k.push_back(s2[k]);
            i += 1;
            k += 1;
            break;

        case UnpairedJL:
            list_j.insert(0, 1, s1[j]);
            list_bp_right.insert(0, 1, '.');
            list_l.insert(0, 1, s2[l]);
            j -= 1;
            l -= 1;
            break;

        case PairedGapS1:
            list_i.push_back(s1[i]);
            list_j.insert(0, 1, s1[j]);
            list_bp_left.push_back('(');
            list_bp_right.insert(0, 1, ')');
            list_k.push_back('-');
            list_l.insert(0, 1, '-');
            i += 1;
            j -= 1;
            break;

        case PairedGapS2:
            list_k.push_back(s2[k]);
            list_l.insert(0, 1, s2[l]);
            list_bp_left.push_back('(');
            list_bp_right.insert(0, 1, ')');
            list_i.push_back('-');
            list_j.insert(0, 1, '-');
            k += 1;
            l -= 1;
            break;

        case Paired:
            list_i.push_back(s1[i]);
            list_j.insert(0, 1, s1[j]);
            list_bp_left.push_back('(');
            list_bp_right.insert(0, 1, ')');
            list_k.push_back(s2[k]);
            list_l.insert(0, 1, s2[l]);
            i += 1;
            j -= 1;
            k += 1;
            l -= 1;
            break;
    }
    return dp_matrix->get_pos(i, j, k, l);
}

void Backtrace::calculate_mb_position(float score)
{
    for (m = i + 1; m < j; ++m)
    {
        for (n = k + 1; n < l; ++n)
        {
            dp_matrix_cell mb_right = dp_matrix->get_pos(m + 1, j, n + 1, l);
            if (mb_right.parent != Paired)
                continue;

            if (dp_matrix->get_pos(i, m, k, n).score + mb_right.score == score)
                return;
        } //n
    } //m
}

void Backtrace::do_backtrace_mb(int i, int j, int k, int l)
{
    std::string temp_s1, temp_structure, temp_s2;
    Backtrace *mb_bc = new Backtrace(dp_matrix, i, j, k, l, s1, s2);
    mb_bc->run();
    mb_bc->print(temp_s1, temp_structure, temp_s2);

    list_i.append(temp_s1);
    list_bp_left.append(temp_structure);
    list_k.append(temp_s2);
    delete mb_bc;
}

void Backtrace::run()
{
    dp_matrix_cell c = dp_matrix->get_pos(i, j, k, l);
    while (dp_matrix->check_border(i, j, k, l) && c.parent != Multibranch)
    {
        //printf("%f (%s) - %d %d %d %d\n", c.score, parent_str[(int)c.parent], i, k, j, l);
        c = get_parent(c);
    }
    if (c.parent == Multibranch)
    {
        calculate_mb_position(c.score);
        //printf("(%d %d %d %d) = (%d %d %d %d) + (%d %d %d %d)\n", i, j, k, l, i, m, k, n, m + 1, j, n + 1, l);

        do_backtrace_mb(i, m, k, n);
        do_backtrace_mb(m + 1, j, n + 1, l);
        return;
    }
    //printf("Fim: %f (%s) - %d %d %d %d\n", c.score, parent_str[(int)c.parent], i, k, j, l);
}

void Backtrace::print(std::string &alignment_s1, std::string &alignment_structure, std::string &alignment_s2)
{
    alignment_s1.append(list_i);
    alignment_s1.append(list_j);
    alignment_structure.append(list_bp_left);
    alignment_structure.append(list_bp_right);
    alignment_s2.append(list_k);
    alignment_s2.append(list_l);
}
