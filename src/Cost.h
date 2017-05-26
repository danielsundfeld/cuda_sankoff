#ifndef _COST_H
#define _COST_H
#ifdef __CUDACC__
#define CUDAFLAGS __device__ __host__
#else
#define CUDAFLAGS
#endif

class Cost {
    public:
        static const float gap = -3;
        static const float unpaired = 0.05;

        CUDAFLAGS static int base_score(const char &a, const char &b); //'alpha' score
        CUDAFLAGS static float match_score(const char &a, const char &b); // 'beta' score
        CUDAFLAGS static int compensation_score(const char &a1, const char &a2, const char &b1, const char &b2); // 'tau' score
};

CUDAFLAGS inline int Cost::base_score(const char &a, const char &b)
{
    if ((a == 'C' && b == 'G') ||
        (a == 'G' && b == 'C') ||
        (a == 'A' && b == 'U') ||
        (a == 'U' && b == 'A') ||
        (a == 'G' && b == 'U') ||
        (a == 'U' && b == 'G')
       )
        return 1;
    return 0;
}

//Only for unpaired
CUDAFLAGS inline float Cost::match_score(const char &a, const char &b)
{
    if (a == b)
        return Cost::unpaired;
    return 0;
}

CUDAFLAGS inline int Cost::compensation_score(const char &a1, const char &a2, const char &b1, const char &b2)
{
    // Give a score if they have the same structure
    if (base_score(a1, a2) == base_score(b1, b2))
        return 0; //TODO: FIXME 2 * match; //4 nucleotides, 2 matches...
    return 0; // Otherwise, no donuts for you
}
#endif
