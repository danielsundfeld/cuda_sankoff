#ifndef _COST_H
#define _COST_H
#ifdef __CUDACC__
#define CUDAFLAGS __device__ __host__
#else
#define CUDAFLAGS
#endif

#define _RIBOSUM_8560 1

class Cost {
    public:
        static const float gap = -13;
        static const float unpaired = 5;
        static const float match = 0.05;

        CUDAFLAGS static int base_score(const char &a, const char &b); //'alpha' score
        CUDAFLAGS static float match_score(const char &a, const char &b); // 'beta' score
        CUDAFLAGS static float compensation_score(const char &a1, const char &a2, const char &b1, const char &b2); // 'tau' score
};

CUDAFLAGS inline int Cost::base_score(const char &a, const char &b)
{
    /* In the BP scores, this function only checks if it is a base pair */
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

#ifdef _RIBOSUM_8560
CUDAFLAGS inline float Cost::match_score(const char &a, const char &b)
{
    int x, y;
    const float ribosum_85_60_substitution[4][4] = {
        //A     C      G      U
        { 2.22, -1.86, -1.46, -1.39 }, //A
        {-1.86,  1.16, -2.48, -1.05 }, //C
        {-1.46, -2.48,  1.03, -1.74 }, //G
        {-1.39, -1.05, -1.74,  1.65 }, //U
    };

    if (a == 'A')
        x = 0;
    else if (a == 'C')
        x = 1;
    else if (a == 'G')
        x = 2;
    else if (a == 'U')
        x = 3;
    else
        return -1024;

    if (b == 'A')
        y = 0;
    else if (b == 'C')
        y = 1;
    else if (b == 'G')
        y = 2;
    else if (b == 'U')
        y = 3;
    else
        return -1024;

    return ribosum_85_60_substitution[x][y];
}

CUDAFLAGS inline float Cost::compensation_score(const char &a1, const char &a2, const char &b1, const char &b2)
{
    const float ribosum_85_60_compensation[16][16] = {
        //  AA       AC,     AG,     AU,     CA,     CC,     CG,     CU,     GA,     GC,     GG,     GU,     UA,     UC,     UG,    UU
        { -2.49,  -7.04,  -8.24,  -4.32,  -8.84, -14.37,  -4.68, -12.64,  -6.86,  -5.03,  -8.39,  -5.84,  -4.01, -11.32,  -6.16,  -9.05}, //AA
        { -7.04,  -2.11,  -8.89,  -2.04,  -9.37,  -9.08,  -5.86, -10.45,  -9.73,  -3.81, -11.05, -4.72 ,  -5.33,  -8.67,  -6.93,  -7.83}, //AC
        { -8.24,  -8.89,  -0.80,  -5.13, -10.41, -14.53,  -4.57, -10.14,  -8.61,  -5.77, -5.38 , -6.60 ,  -5.43,  -8.87, -5.94 , -11.07}, //AG
        { -4.32,  -2.04,  -5.13,   4.49,  -5.56,  -6.71,   1.67,  -5.17,  -5.33,   2.70, -5.61 ,  0.59 ,   1.61, -4.81 ,  -0.51,  -2.98}, //AU
        { -8.84,  -9.37, -10.41,  -5.56,  -5.13, -10.45,  -3.57,  -8.49,  -7.98,  -5.95, -11.36, -7.93 ,  -2.42,  -7.08,  -5.63,  -8.39}, //CA
        {-14.37,  -9.08, -14.53,  -6.71, -10.45,  -3.59,  -5.71,  -5.77, -12.43,  -3.70, -12.58,  -7.88,  -6.88,  -7.40,  -8.41,  -5.41}, //CC
        { -4.68,  -5.86,  -4.57,   1.67,  -3.57,  -5.71,   5.36,  -4.96,  -6.00,   2.11,  -4.66,  -0.27,   2.75,  -4.91,  1.32 ,  -3.67}, //CG
        {-12.64, -10.45, -10.14, -5.17 ,  -8.49,  -5.77,  -4.96,  -2.28,  -7.71,  -5.84, -13.69,  -5.61,  -4.72,  -3.83,  -7.36,  -5.21}, //CU
        {-6.86 ,  -9.73,  -8.61,  -5.33,  -7.98,-12.43 ,  -6.00,  -7.71,  -1.05,  -4.88,  -8.67, -6.10 ,  -5.85,  -6.63,  -7.55, -11.54}, //GA
        { -5.03,  -3.81,  -5.77,   2.70,  -5.95,  -3.70,   2.11,  -5.84,  -4.88,   5.62,  -4.13,  1.21 ,  1.60 ,  -4.49,  -0.08,  -3.90}, //GC
        { -8.39, -11.05,  -5.38,  -5.61, -11.36, -12.58,  -4.66, -13.69,  -8.67,  -4.13,  -1.98,  -5.77,  -5.75, -12.01,  -4.27, -10.79}, //GG
        { -5.84,  -4.72,  -6.60,   0.59,  -7.93,  -7.88, -0.27 ,  -5.61,  -6.10,   1.21,  -5.77,   3.47,  -0.57,  -5.30,  -2.09,  -4.45}, //GU
        { -4.01,  -5.33,  -5.43,   1.61,  -2.42,  -6.88,   2.75,  -4.72,  -5.85,   1.60,  -5.75,  -0.57,   4.97,  -2.98,   1.14,  -3.39}, //UA
        {-11.32,  -8.67,  -8.87,  -4.81,  -7.08,  -7.40,  -4.91,  -3.83,  -6.63,  -4.49, -12.01,  -5.30,  -2.98,  -3.21,  -4.76,  -5.97}, //UC
        { -6.16,  -6.93,  -5.94,  -0.51,  -5.63,  -8.41,   1.32,  -7.36,  -7.55,  -0.08,  -4.27,  -2.09,   1.14,  -4.76,   3.36,  -4.28}, //UG
        { -9.05,  -7.83, -11.07,  -2.98,  -8.39,  -5.41,  -3.67,  -5.21, -11.54,  -3.90, -10.79,  -4.45,  -3.39,  -5.97,  -4.28,  -0.02}, //UU
    };
    int x, y;
    if (a1 == 'A' && a2 == 'A')
        x = 0;
    else if (a1 == 'A' && a2 == 'C')
        x = 1;
    else if (a1 == 'A' && a2 == 'G')
        x = 2;
    else if (a1 == 'A' && a2 == 'U')
        x = 3;
    else if (a1 == 'C' && a2 == 'A')
        x = 4;
    else if (a1 == 'C' && a2 == 'C')
        x = 5;
    else if (a1 == 'C' && a2 == 'G')
        x = 6;
    else if (a1 == 'C' && a2 == 'U')
        x = 7;
    else if (a1 == 'G' && a2 == 'A')
        x = 8;
    else if (a1 == 'G' && a2 == 'C')
        x = 9;
    else if (a1 == 'G' && a2 == 'G')
        x = 10;
    else if (a1 == 'G' && a2 == 'U')
        x = 11;
    else if (a1 == 'U' && a2 == 'A')
        x = 12;
    else if (a1 == 'U' && a2 == 'C')
        x = 13;
    else if (a1 == 'U' && a2 == 'G')
        x = 14;
    else if (a1 == 'U' && a2 == 'U')
        x = 15;
    else
        return -1024;

    if (b1 == 'A' && b2 == 'A')
        y = 0;
    else if (b1 == 'A' && b2 == 'C')
        y = 1;
    else if (b1 == 'A' && b2 == 'G')
        y = 2;
    else if (b1 == 'A' && b2 == 'U')
        y = 3;
    else if (b1 == 'C' && b2 == 'A')
        y = 4;
    else if (b1 == 'C' && b2 == 'C')
        y = 5;
    else if (b1 == 'C' && b2 == 'G')
        y = 6;
    else if (b1 == 'C' && b2 == 'U')
        y = 7;
    else if (b1 == 'G' && b2 == 'A')
        y = 8;
    else if (b1 == 'G' && b2 == 'C')
        y = 9;
    else if (b1 == 'G' && b2 == 'G')
        y = 10;
    else if (b1 == 'G' && b2 == 'U')
        y = 11;
    else if (b1 == 'U' && b2 == 'A')
        y = 12;
    else if (b1 == 'U' && b2 == 'C')
        y = 13;
    else if (b1 == 'U' && b2 == 'G')
        y = 14;
    else if (b1 == 'U' && b2 == 'U')
        y = 15;
    else
        return -1024;

    return ribosum_85_60_compensation[x][y];
}
#else
CUDAFLAGS inline float Cost::match_score(const char &a, const char &b)
{
    if (a == b)
        return Cost::match;
    return 0;
}

CUDAFLAGS inline float Cost::compensation_score(const char &a1, const char &a2, const char &b1, const char &b2)
{
    // Give a score if they have the same structure
    if (base_score(a1, a2) == base_score(b1, b2))
        return Cost::match;
    return 0; // Otherwise, no donuts for you
}
#endif
#endif
