#include "Sankoff.h"

#include <iostream>
#include <stdio.h>

#include "bp_probs.h"
#include "Cost.h"
#include "TimeCounter.h"

Sankoff::Sankoff(const std::string &seq1, const std::string &seq2)
: dp_matrix(seq1.length(), seq2.length())
{
    s1 = seq1;
    s1_l = (int) s1.length(); //TODO: throw exception if negative
    s2 = seq2;
    s2_l = (int) s2.length();
    bp1 = new bp_prob();
    bp2 = new bp_prob();

    get_bp_prob(seq1, bp1);
    get_bp_prob(seq2, bp2);
}

Sankoff::~Sankoff()
{
    delete bp1;
    delete bp2;
}

void Sankoff::print_orig(int i, int j, int k, int l) const
{
    std::cout << i << " " << k << " " << j << " "  << l << ":";
}

void Sankoff::print_index(int i, int j, int k, int l) const
{
    if (dp_matrix.check_border(i,j,k,l))
    {
        std::cout << "\t" << i << " " << k << " " << j << " "  << l << " " <<
        dp_matrix.get_pos(i, j, k, l).score;
    }
    std::cout << "\n";
}

// To visualize the data dependecy...
void Sankoff::print_score_dep(int i, int j, int k, int l) const
{
    return;
    print_orig(i, j, k, l);
    std::cout << " (GapI) ";
    print_index(i + 1, j, k, l);
    print_orig(i, j, k, l);
    std::cout << " (GapK) ";
    print_index(i, j, k + 1, l);
    print_orig(i, j, k, l);
    std::cout << " (GapJ) ";
    print_index(i, j - 1, k, l);
    print_orig(i, j, k, l);
    std::cout << " (GapL) ";
    print_index(i, j, k, l - 1);

    print_orig(i, j, k, l);
    std::cout << " (UnpairedIK) ";
    print_index(i + 1, j, k + 1, l);
    print_orig(i, j, k, l);
    std::cout << " (UnpairedJL) ";
    print_index(i, j - 1, k, l - 1);

    print_orig(i, j, k, l);
    std::cout << " (Paired GapS1) ";
    print_index(i + 1, j - 1, k, l);
    print_orig(i, j, k, l);
    std::cout << " (Paired GapS2) ";
    print_index(i, j, k + 1, l - 1);

    print_orig(i, j, k, l);
    std::cout << " (Paired) ";
    print_index(i + 1, j - 1 , k + 1, l - 1);
}

// To visualize the data dependecy...
void Sankoff::print_mb_dep(int i, int j, int k, int l, int m, int n) const
{
    return;
    std::cout << i << " " << j << " " << k << " " << l << ":\t"
        << i << " " << m << " " << k << " " << n
        << " + "
        << m + 1 << " " << j << " " << n + 1 << " " << l << "\n";
}

void Sankoff::max(dp_matrix_cell &score1, dp_matrix_cell score2, int parent)
{
    if (score2.score > score1.score || (score1.parent == NullParent))
    {
        score1.score = score2.score;
        score1.parent = parent;
    }
}

void Sankoff::calculate_pos(dp_matrix_cell &score1, int i, int j, int k, int l, float extra_score, int parent)
{
    if (dp_matrix.check_border(i, j, k, l) == false)
        return;

    dp_matrix_cell score2(dp_matrix.get_pos(i, j, k, l));
    score2.score += extra_score;
    max(score1, score2, parent);
}

void Sankoff::calculate_pos_mb(dp_matrix_cell &score1, int i, int j, int k, int l, int m, int n)
{
    dp_matrix_cell mb_right(dp_matrix.get_pos(m + 1, j, n + 1, l));
    if (mb_right.parent != Paired)
        return;

    dp_matrix_cell mb_left(dp_matrix.get_pos(i, m, k, n));
    dp_matrix_cell mb;
    mb.score = mb_left.score + mb_right.score;
    mb.parent = Multibranch;

    max(score1, mb, Multibranch);
}

void Sankoff::expand_pos(const int &i, const int &j, const int &k, const int &l)
{
    dp_matrix_cell score = dp_matrix_cell();

    if (dp_matrix.check_border(i, j, k, l) == false)
        return;

    print_score_dep(i, j, k, l);
    float s1_score = bp1->m[i+1][j+1];
    float s2_score = bp2->m[k+1][l+1];
    if (s1_score == 0 && s2_score == 0)
    {
        s1_score = -1024;
        s2_score = -1024;
    }

    /*
     * Explanations of this recursion functions can be see at:
     *
     * - Havgaard, et al. "Fast Pairwise Structural RNA Alignments by
     * Pruning of the Dynamical Programming Matrix".
     * - Havgaard, et al. "Pairwise local structural alignment of RNA sequences
     * with sequence similarity less than 40%".
     * - Torarinsson, et al. "Multiple structural alignment and clustering of RNA sequences
     * - Ziv-Ukelson, et al. "A faster algorithm for RNA co-folding"
    */
    calculate_pos(score, i + 1, j, k, l, Cost::gap, GapI);
    calculate_pos(score, i, j, k + 1, l, Cost::gap, GapK);
    calculate_pos(score, i, j - 1, k, l, Cost::gap, GapJ);
    calculate_pos(score, i, j, k, l - 1, Cost::gap, GapL);
    calculate_pos(score, i + 1, j, k + 1, l, Cost::match_score(s1[i], s2[k]), UnpairedIK);
    calculate_pos(score, i, j - 1, k, l - 1, Cost::match_score(s1[j], s2[l]), UnpairedJL);
    calculate_pos(score, i + 1, j - 1, k, l, s1_score + Cost::gap * 2, PairedGapS1);
    calculate_pos(score, i, j, k + 1, l - 1, s2_score + Cost::gap * 2, PairedGapS2);
    calculate_pos(score, i + 1, j - 1, k + 1, l - 1, s1_score + s2_score +
            Cost::compensation_score(s1[i], s1[j], s2[k], s2[l]), Paired);

    for (int m = i + 1; m < j; ++m)
    {
        for (int n = k + 1; n < l; ++n)
        {
            print_mb_dep(i, j, k, l, m, n);
            calculate_pos_mb(score, i, j, k, l, m, n);
        } //n
    } //m

    //printf("%f (%s) %d %d %d %d\n", score.score, parent_str[(int)score.parent], i, k, j, l);
    dp_matrix.put_pos(i, j, k, l, score);
}

void Sankoff::expand_inner_matrix(const int &i, const int &k)
{
    for (int j = i; j < s1_l; ++j)
    {
        for (int l = k; l < s2_l; ++l)
        {
            expand_pos(i, j, k, l);
        } //l
    } //j
}

void Sankoff::expand_inner_matrix_diag(const int &i, const int &k)
{
    if (i < 0 || i >= s1_l || k < 0 || k >= s2_l)
        return;

    for (int inner_diag = i; inner_diag < s1_l; ++inner_diag)
    {
#pragma omp parallel for schedule(dynamic,1)
        for (int l = k; l < s2_l; ++l)
        {
            int j = inner_diag - (l - k);
            expand_pos(i, j, k, l);
        }
    }
    for (int inner_diag = k + 1; inner_diag < s2_l; ++inner_diag)
    {
#pragma omp parallel for schedule(dynamic,1)
        for (int l = inner_diag; l < s1_l; ++l)
        {
            int j = s1_l - 1 - (l - inner_diag);
            expand_pos(i, j, k, l);
        }
    }
}

int Sankoff::sankoff()
{
    std::cout << "Sankoff:"
        << "\nseq1: (" << s1_l << ")\t" << s1
        << "\nseq2: (" << s2_l << ")\t" << s2
        << "\n";

    for (int i = s1_l - 1; i >= 0; --i)
    {
        for (int k = s2_l - 1; k >= 0; --k)
        {
            expand_inner_matrix(i, k);
        } //k
    } //i
    std::cout << dp_matrix.get_pos(0, s1_l - 1, 0, s2_l - 1).score << std::endl;
    return 0;
}

int Sankoff::diag_sankoff()
{
    std::cout << "Sankoff:"
        << "\nseq1: (" << s1_l << ")\t" << s1
        << "\nseq2: (" << s2_l << ")\t" << s2
        << "\n";

    for (int outer_diag = 0; outer_diag <= s2_l - 1; ++outer_diag)
    {
#pragma omp parallel for schedule(dynamic,1)
        for (int k = s2_l - 1 - outer_diag; k <= s2_l - 1; ++k)
        {
            int i = s2_l - 1 - outer_diag + s1_l - 1 - k;
            expand_inner_matrix_diag(i, k);
        } //i
    } //outer_diag
    for (int outer_diag = s1_l - 2; outer_diag >= 0 ; --outer_diag)
    {
#pragma omp parallel for schedule(dynamic,1)
        for (int k = 0; k <= s2_l - 1; ++k)
        {
            int i = outer_diag - k;
            if (i >= 0)
                expand_inner_matrix_diag(i, k);
        }
    } //outer_diag
    std::cout << "Score: " << dp_matrix.get_pos(0, s1_l - 1, 0, s2_l - 1).score << std::endl;
    //dp_matrix.backtrace("1234567890123456789", "1234567890123456789");
    TimeCounter tb("Backtrace total time");
    dp_matrix.backtrace(s1, s2);
    return 0;
}
