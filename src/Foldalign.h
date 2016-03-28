#ifndef _FOLDALIGN_H
#define _FOLDALIGN_H
#include <string>

#include "Sankoff.h"

class Foldalign : public Sankoff {
    public:
        Foldalign(const std::string &seq1, const std::string &seq2, const int lmbd, const int dlt);
        int fold_align(); //Run with the Foldalign heuristic

    private:
        int lambda;
        int delta;
};
#endif
