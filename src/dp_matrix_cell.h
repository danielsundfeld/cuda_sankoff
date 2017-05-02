#ifndef _DP_MATRIX_CELL_H
#define _DP_MATRIX_CELL_H
struct dp_matrix_cell {
    float score;
    char parent;
};

enum Parent {
    NullParent = 0,
    FirstCell,
    GapI,
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

const char parent_str[][LastParent] = {
    "NullParent",
    "FirstCell",
    "GapI",
    "GapK",
    "GapJ",
    "GapL",
    "UnpairedIK",
    "UnpairedJL",
    "PairedGapS1",
    "PairedGapS2",
    "Paired",
    "Multibranch",
    "LastParent"
};
#endif
