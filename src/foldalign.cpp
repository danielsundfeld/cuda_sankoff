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

int Foldalign::fold_align()
{
    std::cout << "Foldalign:"
        << "\n\u03B4 = " << m_delta
        << "\n\u03BB = " << m_lambda
        << "\nseq1:\t" << m_seq1
        << "\nseq2:\t" << m_seq2
        << "\n";

    return 0;
}
