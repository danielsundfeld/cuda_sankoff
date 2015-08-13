#include "Cost.h"

int Cost::base_score(const char &a, const char &b)
{
    if ((a == 'C' && b == 'G') ||
        (a == 'G' && b == 'C') ||
        (a == 'A' && b == 'U') ||
        (a == 'U' && b == 'A')
       )
        return base_paired;
    return unpaired;
}

int Cost::match_score(const char &a, const char &b)
{
    if (a == b)
        return match;
    return mismatch;
}

int Cost::compensation_score(const char &a1, const char &a2, const char &b1, const char &b2)
{
    if (match_score(a1, b1) && match_score(a2, b2))
        return unpaired; // If equal, compensation score is zero
    
    if (base_score(a1, a2) && base_score(b1, b2))
        return match * 2; // They are different, but it is a basepair! compensate

    return unpaired;
}
