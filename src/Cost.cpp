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
    /*
     * The score function is alpha + beta + tao
     * alpha (a): base_score
     * beta (b): match_score
     * tao (t): compensation_score
     */

    // If equal, a = match, b = base_pair, t = zero
    if (match_score(a1, b1) == match && match_score(a2, b2) == match)
        return 0;
    
    /* They are different but they are basepaired and have the same structure
       a = mismatch, b = base_pair, t = match - mismatch (compensate) */
    if (base_score(a1, a2) == base_paired && base_score(b1, b2) == base_paired)
        return match - mismatch;

    return 0; // Otherwise, no donuts for you
}
