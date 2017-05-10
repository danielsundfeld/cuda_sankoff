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

        case Multibranch:
            std::cout << "MULTIBRANCH\n";
            //exit(1);
            break;
    }
    return dp_matrix->get_pos(i, j, k, l);
}

void Backtrace::add_last(const dp_matrix_cell c)
{
    if (c.parent != NullParent)
        return;

    list_i.push_back(s1[i]);
    list_k.push_back(s2[k]);
    list_bp_left.push_back('.');
}

void print_list(const std::string &list, std::string &st)
{
    st.append(list);
}

void Backtrace::run()
{
    dp_matrix_cell c = dp_matrix->get_pos(i, j, k, l);
    while (dp_matrix->check_border(i, j, k, l) && c.parent != NullParent)
    {
        if (c.parent == Multibranch)
            return;
        printf("%f (%s) - %d %d %d %d\n", c.score, parent_str[(int)c.parent], i, k, j, l);
        c = get_parent(c);
    }
    add_last(c);
    printf("Fim: %f (%s) - %d %d %d %d\n", c.score, parent_str[(int)c.parent], i, k, j, l);
}

void Backtrace::print(std::string &alignment_s1, std::string &alignment_structure, std::string &alignment_s2)
{
    print_list(list_i, alignment_s1);
    print_list(list_j, alignment_s1);
    print_list(list_bp_left, alignment_structure);
    print_list(list_bp_right, alignment_structure);
    print_list(list_k, alignment_s2);
    print_list(list_l, alignment_s2);
}
