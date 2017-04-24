#ifndef _DP_MATRIX_CELL_H
#define _DP_MATRIX_CELL_H
struct dp_matrix_cell {
    float score;
    char parent;
};

enum Parent {
    GapI = 0,
    GapK,
    GapJ,
    GapL,
    UnpairedIK,
    UnpairedJL,
    PairedGapS1,
    PairedGapS2,
    Paired,
    Multibranch,
    LastParent
};
#endif
