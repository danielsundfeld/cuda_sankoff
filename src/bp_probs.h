#ifndef _BP_PROBS
#define _BP_PROBS
#include <string>

#include "DPMatrix_GPU.h"

struct bp_prob {
    int size;
    float m[MAX_SEQ_SIZE][MAX_SEQ_SIZE];
};

int get_bp_prob(const std::string &seq, struct bp_prob *bp);
#endif
