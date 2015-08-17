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
    // Give a score if they have the same structure
    if (base_score(a1, a2) == base_score(b1, b2))
        return 2 * match; //4 nucleotides, 2 matches...
    return 0; // Otherwise, no donuts for you
}
