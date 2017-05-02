#ifndef _COST_H
#define _COST_H
#ifdef __CUDACC__
#define CUDAFLAGS __device__ __host__
#else
#define CUDAFLAGS
#endif

class Cost {
    public:
        static const float gap = -0.4;
        static const float mismatch = -0.3;
        static const float match = 0.3;
        static const float base_paired = 0.3;
        static const float unpaired = -0.2;

        CUDAFLAGS static int base_score(const char &a, const char &b); //'alpha' score
        CUDAFLAGS static int match_score(const char &a, const char &b); // 'beta' score
        CUDAFLAGS static int compensation_score(const char &a1, const char &a2, const char &b1, const char &b2); // 'tau' score
};

CUDAFLAGS inline int Cost::base_score(const char &a, const char &b)
{
    if ((a == 'C' && b == 'G') ||
        (a == 'G' && b == 'C') ||
        (a == 'A' && b == 'U') ||
        (a == 'U' && b == 'A')
       )
        return base_paired;
    return unpaired;
}

CUDAFLAGS inline int Cost::match_score(const char &a, const char &b)
{
    if (a == b)
        return match;
    return mismatch;
}

CUDAFLAGS inline int Cost::compensation_score(const char &a1, const char &a2, const char &b1, const char &b2)
{
    // Give a score if they have the same structure
    if (base_score(a1, a2) == base_score(b1, b2))
        return 2 * match; //4 nucleotides, 2 matches...
    return 0; // Otherwise, no donuts for you
}
#endif
