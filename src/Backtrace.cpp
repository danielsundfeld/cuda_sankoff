#include "Backtrace.h"

#include <iostream>
#include <string>
#include <stdio.h>

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

void print_list(const std::list<char> &list)
{
    for (std::list<char>::const_iterator it = list.begin(); it != list.end(); ++it)
        std::cout << *it;
}

void Backtrace::run()
{
    dp_matrix_cell c = dp_matrix->get_pos(i, j, k, l);
    while (dp_matrix->check_border(i, j, k, l) && c.parent != NullParent)
    {
        printf("%f (%s) - %d %d %d %d\n", c.score, parent_str[(int)c.parent], i, k, j, l);
        c = get_parent(c);
    }
    add_last(c);
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
