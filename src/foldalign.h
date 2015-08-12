#ifndef _FOLDALIGN_H
#define _FOLDALIGN_H
#include <string>

class Foldalign {
    public:
        Foldalign(const std::string &s1, const std::string &s2, const int lambda, const int delta);
        int fold_align();
    private:
        std::string m_seq1;
        std::string m_seq2;
        int m_lambda;
        int m_delta;
};
#endif
