#ifndef _COST_H
#define _COST_H
class Cost {
    public:
        enum {
            gap = 0,
            mismatch = 1,
            match = 2,

            base_paired = 2,
            unpaired = 0,
        };
        static int base_score(const char &a, const char &b); //'alpha' score
        static int match_score(const char &a, const char &b); // 'beta' score
        static int compensation_score(const char &a1, const char &a2, const char &b1, const char &b2); // 'tau' score
};
#endif
