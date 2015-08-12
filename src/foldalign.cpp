#include "foldalign.h"

#include <iostream>

Foldalign::Foldalign(std::string s1, std::string s2)
{
    m_seq1 = s1;
    m_seq2 = s2;
}

int Foldalign::fold_align()
{
    std::cout << "Folding:\n\t" << m_seq1 << "\n\t" << m_seq2 << "\n";
    return 0;
}
