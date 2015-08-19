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

void Foldalign::print_index(int i, int j, int k, int l) const
{
    std::cout << "\t" << i << " " << j << " " << k << " "  << l << "\n";
}

// To visualize the data dependecy...
void Foldalign::print_score_dep(int i, int j, int k, int l) const
{
    print_orig(i, j, k, l);
    print_index(i + 1, j, k, l);
    print_orig(i, j, k, l);
    print_index(i, j, k + 1, l);
    print_orig(i, j, k, l);
    print_index(i, j - 1, k, l);
    print_orig(i, j, k, l);
    print_index(i, j, k, l - 1);

    print_orig(i, j, k, l);
    print_index(i + 1, j, k + 1, l);
    print_orig(i, j, k, l);
    print_index(i, j - 1, k, l - 1);

    print_orig(i, j, k, l);
    print_index(i + 1, j - 1, k, l);
    print_orig(i, j, k, l);
    print_index(i, j, k + 1, l - 1);

    print_orig(i, j, k, l);
    print_index(i + 1, j - 1 , k + 1, l - 1);
}

// To visualize the data dependecy...
void Foldalign::print_mb_dep(int i, int j, int k, int l, int m, int n) const
{
    std::cout << i << " " << j << " " << k << " " << l << ":\t"
        << i << " " << m << " " << k << " " << n
        << " + "
        << m + 1 << " " << j << " " << n + 1 << " " << l << "\n";
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

bool Foldalign::out_of_border_lambda(const int j) const
{
    if (j >= m_seq1_l)
        return true;
    return false;
}

bool Foldalign::out_of_border_delta(const int l) const
{
    if (l < 0 || l >= m_seq2_l)
        return true;
    return false;
}

int Foldalign::sankoff()
{
    std::cout << "Sankoff:"
        << "\nseq1:\t" << m_seq1
        << "\nseq2:\t" << m_seq2
        << "\n";

    for (int i = m_seq1_l - 1; i >= 0; --i)
    {
        for (int k = m_seq2_l - 1; k >= 0; --k)
        {
            for (int j = i + 1; j < m_seq1_l; ++j)
            {
                for (int l = k + 1; l < m_seq2_l; ++l)
                {
                    int score = 0;

                    print_score_dep(i, j, k, l);

                    score = std::max(score, dp_matrix[index(i + 1, j, k, l)] + Cost::gap);
                    score = std::max(score, dp_matrix[index(i, j, k + 1, l)] + Cost::gap);
                    score = std::max(score, dp_matrix[index(i, j - 1, k, l)] + Cost::gap);
                    score = std::max(score, dp_matrix[index(i, j, k, l - 1)] + Cost::gap);
                    score = std::max(score, dp_matrix[index(i + 1, j, k + 1, l)] + Cost::match_score(m_seq1[i], m_seq2[k]));
                    score = std::max(score, dp_matrix[index(i, j - 1, k, l - 1)] + Cost::match_score(m_seq1[j], m_seq2[l]));
                    score = std::max(score, dp_matrix[index(i + 1, j - 1, k, l)] + Cost::base_score(m_seq1[i], m_seq1[j]) + Cost::gap * 2);
                    score = std::max(score, dp_matrix[index(i, j, k + 1, l - 1)] + Cost::base_score(m_seq2[k], m_seq2[l]) + Cost::gap * 2);
                    score = std::max(score, dp_matrix[index(i + 1, j - 1, k + 1, l - 1)] +
                                                        Cost::base_score(m_seq1[i], m_seq1[j]) + Cost::base_score(m_seq2[k], m_seq2[l]) +
                                                        Cost::compensation_score(m_seq1[i], m_seq1[j], m_seq2[k], m_seq2[l]));

                    for (int m = i + 1; m < j; ++m)
                    {
                        for (int n = k + 1; n < l; ++n)
                        {
                            print_mb_dep(i, j, k, l, m, n);
                            score = std::max(score, dp_matrix[index(i, m, k, n)] + dp_matrix[index(m + 1, j, n + 1, l)]);
                        } //n
                    } //m
                    if (score > 0)
                        dp_matrix[index(i, j, k, l)] = score;
                } //l
            } //j
        } //k
    } //i
    std::cout << dp_matrix[index(0, m_seq1_l - 1, 0, m_seq2_l - 1)] << std::endl;
    return 0;
}

int Foldalign::fold_align()
{
    std::cout << "Foldalign:"
        << "\n\u03B4 = " << m_delta
        << "\n\u03BB = " << m_lambda
        << "\nseq1:\t" << m_seq1
        << "\nseq2:\t" << m_seq2
        << "\n";

    for (int i = m_seq1_l - 1; i >= 0; --i)
    {
        for (int k = m_seq2_l - 1; k >= 0; --k)
        {
            int j_end = i + 1 + m_lambda;
            if (j_end > m_seq1_l)
                j_end = m_seq1_l;

            for (int j = i + 1; j < j_end; ++j)
            {
                int l_begin = j - m_delta;
                if (l_begin < k + 1)
                    l_begin = k + 1;
                int l_end = j + m_delta;
                if (l_end > m_seq2_l)
                    l_end = m_seq2_l;

                for (int l = l_begin; l < l_end; ++l)
                {
                    int score = 0;

                    print_score_dep(i, j, k, l);

                    score = std::max(score, dp_matrix[index(i + 1, j, k, l)] + Cost::gap);
                    score = std::max(score, dp_matrix[index(i, j, k + 1, l)] + Cost::gap);
                    score = std::max(score, dp_matrix[index(i, j - 1, k, l)] + Cost::gap);
                    score = std::max(score, dp_matrix[index(i, j, k, l - 1)] + Cost::gap);
                    score = std::max(score, dp_matrix[index(i + 1, j, k + 1, l)] + Cost::match_score(m_seq1[i], m_seq2[k]));
                    score = std::max(score, dp_matrix[index(i, j - 1, k, l - 1)] + Cost::match_score(m_seq1[j], m_seq2[l]));
                    score = std::max(score, dp_matrix[index(i + 1, j - 1, k, l)] + Cost::base_score(m_seq1[i], m_seq1[j]) + Cost::gap * 2);
                    score = std::max(score, dp_matrix[index(i, j, k + 1, l - 1)] + Cost::base_score(m_seq2[k], m_seq2[l]) + Cost::gap * 2);
                    score = std::max(score, dp_matrix[index(i + 1, j - 1, k + 1, l - 1)] +
                                                        Cost::base_score(m_seq1[i], m_seq1[j]) + Cost::base_score(m_seq2[k], m_seq2[l]) +
                                                        Cost::compensation_score(m_seq1[i], m_seq1[j], m_seq2[k], m_seq2[l]));

                    int m_end = i + 1 + m_lambda;
                    if (m_end > j)
                        m_end = j;
                    for (int m = i + 1; m < m_end; ++m)
                    {
                        int n_begin = m - m_delta;
                        if (n_begin < k + 1)
                            n_begin = k + 1;
                        int n_end = m + m_delta;
                        if (n_end > l)
                            n_end = l;

                        for (int n = n_begin; n < n_end; ++n)
                        {
                            print_mb_dep(i, j, k, l, m, n);
                            score = std::max(score, dp_matrix[index(i, m, k, n)] + dp_matrix[index(m + 1, j, n + 1, l)]);
                        } //n
                    } //m
                    if (score > 0)
                        dp_matrix[index(i, j, k, l)] = score;
                } //l
            } //j
        } //k
    } //i
    std::cout << dp_matrix[index(0, m_seq1_l - 1, 0, m_seq2_l - 1)] << std::endl;
    return 0;
}
