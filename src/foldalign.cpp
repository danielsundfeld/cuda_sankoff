#include "foldalign.h"

#include <iostream>

Foldalign::Foldalign(const std::string &s1, const std::string &s2,
                     const int lambda, const int delta)
{
    m_seq1 = s1;
    m_seq2 = s2;
    m_lambda = lambda;
    m_delta = delta;
}

bool Foldalign::outOfBorder_j(const int j)
{
    if ((unsigned int)j >= m_seq1.length())
        return true;
    return false;
}

bool Foldalign::outOfBorder_l(const int l)
{
    if (l < 0 || (unsigned int)l >= m_seq2.length())
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

    for (int i = m_seq1.length() - 1; i >= 0; --i)
    {
        for (int k = m_seq2.length() - 1; k >= 0; --k)
        {
            for (unsigned int j = 0; j < m_seq1.length(); ++j) //TODO: lambda
            {
                if (outOfBorder_j(j))
                    continue;
                for (unsigned int l = 0; l < m_seq2.length(); ++l) //TODO: delta
                {
                    if (outOfBorder_l(l))
                        continue;
                    //TODO: score function
                    for (int m = j - 1; m < i + 1; ++m) //TODO: lambda
                    {
                        if (outOfBorder_j(j))
                            continue;
                        for (int n = l - 1; n < k + 1; ++n) //TODO: delta
                        {
                            if (outOfBorder_l(l))
                                continue;
                            //TODO: multibranch function
                        }
                    }
                    //std::cout << i << " " << k << " " << j << " " << l << "\n";
                } //l
            } //j
        } //k
    } // i
    return 0;
}
