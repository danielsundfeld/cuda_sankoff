#ifndef _FOLDALIGN_H
#define _FOLDALIGN_H
#include <string>

class Foldalign {
    public:
        Foldalign(std::string s1, std::string s2);
        int fold_align();
    private:
        std::string m_seq1;
        std::string m_seq2;
};
#endif
