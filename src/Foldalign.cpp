#include "Foldalign.h"

#include <iostream>

#include "Cost.h"

Foldalign::Foldalign(const std::string &s1, const std::string &s2,
                     const int lambda, const int delta)
{
    m_seq1 = s1;
    m_seq1_l = (int) m_seq1.length(); //TODO: throw exception if negative
    m_seq2 = s2;
    m_seq2_l = (int) m_seq2.length();
    m_lambda = lambda;
    m_delta = delta;
}

void Foldalign::print_orig(int i, int j, int k, int l) const
{
    std::cout << i << " " << j << " " << k << " "  << l << ":";
}

void Foldalign::print_coord(int i, int j, int k, int l) const
{
    std::cout << "\t" << i << " " << j << " " << k << " "  << l << "\n";
}

int Foldalign::calculate_score(int i, int j, int k, int l)
{
    print_orig(i, j, k, l);
    print_coord(i + 1, j, k, l);
    print_orig(i, j, k, l);
    print_coord(i, j, k + 1, l);
    print_orig(i, j, k, l);
    print_coord(i, j - 1, k, l);
    print_orig(i, j, k, l);
    print_coord(i, j, k, l - 1);

    print_orig(i, j, k, l);
    print_coord(i + 1, j, k + 1, l);
    print_orig(i, j, k, l);
    print_coord(i, j - 1, k, l - 1);

    print_orig(i, j, k, l);
    print_coord(i + 1, j - 1, k, l);
    print_orig(i, j, k, l);
    print_coord(i, j, k + 1, l - 1);

    print_orig(i, j, k, l);
    print_coord(i + 1, j - 1 , k + 1, l - 1);
    return 0;
}

int Foldalign::calculate_mb(int i, int j, int k, int l, int m, int n)
{
    std::cout << i << " " << j << " " << k << " " << l << ":\t"
        << i << " " << m << " " << k << " " << n
        << " + "
        << m + 1 << " " << j << " " << n + 1 << " " << l << "\n";
    return -1;
}

bool Foldalign::out_of_border(const int i, const int j, const int k, const int l) const
{
    if (i < 0 || i >= m_seq1_l)
        return true;
    if (j < 0 || j >= m_seq1_l)
        return true;
    if (k < 0 || k >= m_seq2_l)
        return true;
    if (l < 0 || l >= m_seq2_l)
        return true;
    return false;
}

bool Foldalign::out_of_border_lambda(const int j)
{
    if (j >= m_seq1_l)
        return true;
    return false;
}

bool Foldalign::out_of_border_delta(const int l)
{
    if (l < 0 || l >= m_seq2_l)
        return true;
    return false;
}

int Foldalign::fold_align()
{
    std::cout << "Foldalign:"
        << "\n\u03B4 = " << m_delta
        << "\n\u03BB = " << m_lambda
        << "\nseq1:\t" << m_seq1
        << "\nseq2:\t" << m_seq2
        << "\n";

    for (int j = 1; j < m_seq1_l; ++j) //TODO: lambda
    {
        for (int l = 1; l < m_seq2_l; ++l) //TODO: delta
        {
            for (int i = j - 1; i >= 0; --i)
            {
                for (int k = l - 1; k >= 0; --k)
                {
                    int score = 0;

                    calculate_score(i, j, k, l);

                    score = std::max(score, dp_matrix[coord(i + 1, j, k, l)] + Cost::gap);
                    score = std::max(score, dp_matrix[coord(i, j, k + 1, l)] + Cost::gap);
                    score = std::max(score, dp_matrix[coord(i, j - 1, k, l)] + Cost::gap);
                    score = std::max(score, dp_matrix[coord(i, j, k, l - 1)] + Cost::gap);
                    score = std::max(score, dp_matrix[coord(i + 1, j, k + 1, l)] + Cost::match_score(m_seq1[i], m_seq2[k]));
                    score = std::max(score, dp_matrix[coord(i, j - 1, k, l - 1)] + Cost::match_score(m_seq1[j], m_seq2[l]));
                    score = std::max(score, dp_matrix[coord(i + 1, j - 1, k, l)] + Cost::base_score(m_seq1[i], m_seq1[j]) + Cost::gap * 2);
                    score = std::max(score, dp_matrix[coord(i, j, k + 1, l - 1)] + Cost::base_score(m_seq2[k], m_seq2[l]) + Cost::gap * 2);
                    score = std::max(score, dp_matrix[coord(i + 1, j - 1, k + 1, l - 1)] +
                                                        Cost::base_score(m_seq1[i], m_seq1[j]) + Cost::base_score(m_seq2[k], m_seq2[l]) +
                                                        Cost::compensation_score(m_seq1[i], m_seq1[j], m_seq2[k], m_seq2[l]));

                    for (int m = j - 1; m >= i + 1; --m) //TODO: lambda
                    {
                        for (int n = l - 1; n >= k + 1; --n) //TODO: delta
                        {
                            calculate_mb(i, j, k, l, m, n);
                            score = std::max(score, dp_matrix[coord(i, m, k, n)] + dp_matrix[coord(m + 1, j, n + 1, l)]);
                        } //n
                    } //m
                    if (score > 0)
                        dp_matrix[coord(i, j, k, l)] = score;
                    std::cout << std::endl;
                } //k
            } //i
        } //l
    } //j
    return 0;
}
