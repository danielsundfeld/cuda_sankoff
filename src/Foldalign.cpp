#include "Foldalign.h"

#include <iostream>

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

int Foldalign::calculate_score(int i, int j, int k, int l)
{
    return -1;
}

int Foldalign::calculate_mb(int i, int j, int k, int l, int m, int n)
{
    return -1;
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

    dp_matrix[coord(m_seq1_l - 1, 0, m_seq2_l - 1, 0)] = 0;
    for (int i = m_seq1_l - 1; i >= 0; --i)
    {
        for (int k = m_seq2_l - 1; k >= 0; --k)
        {
            for (int j = 0; j < m_seq1_l; ++j) //TODO: lambda
            {
                if (out_of_border_lambda(j))
                    continue;

                for (int l = 0; l < m_seq2_l; ++l) //TODO: delta
                {
                    int score;

                    if (out_of_border_delta(l))
                        continue;
                    score = calculate_score(i, j, k, l);
                    for (int m = j - 1; m < i + 1; ++m) //TODO: lambda
                    {
                        if (out_of_border_lambda(j))
                            continue;
                        for (int n = l - 1; n < k + 1; ++n) //TODO: delta
                        {
                            if (out_of_border_delta(l))
                                continue;
                            score = std::max(score, calculate_mb(i, j, k, l, m, n));
                        } // n
                    } // m
                    dp_matrix[coord(i, j, k, l)] = score;
                    //std::cout << i << " " << k << " " << j << " " << l << "\n";
                } //l
            } //j
        } //k
    } // i
    return 0;
}
