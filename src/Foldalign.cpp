#include "Foldalign.h"

#include <iostream>

#include "Cost.h"

Foldalign::Foldalign(const std::string &seq1, const std::string &seq2,
                     const int lmbd, const int dlt)
{
    s1 = seq1;
    s1_l = (int) s1.length(); //TODO: throw exception if negative
    s2 = seq2;
    s2_l = (int) s2.length();
    lambda = lmbd;
    delta = dlt;
}

void Foldalign::print_orig(int i, int j, int k, int l) const
{
    std::cout << i << " " << j << " " << k << " "  << l << ":";
}

void Foldalign::print_index(int i, int j, int k, int l) const
{
    std::cout << "\t" << i << " " << j << " " << k << " "  << l << "\n";
}

// To visualize the data dependecy...
void Foldalign::print_score_dep(int i, int j, int k, int l) const
{
    return;
    print_orig(i, j, k, l);
    print_index(i + 1, j, k, l);
    print_orig(i, j, k, l);
    print_index(i, j, k + 1, l);
    print_orig(i, j, k, l);
    print_index(i, j - 1, k, l);
    print_orig(i, j, k, l);
    print_index(i, j, k, l - 1);

    print_orig(i, j, k, l);
    print_index(i + 1, j, k + 1, l);
    print_orig(i, j, k, l);
    print_index(i, j - 1, k, l - 1);

    print_orig(i, j, k, l);
    print_index(i + 1, j - 1, k, l);
    print_orig(i, j, k, l);
    print_index(i, j, k + 1, l - 1);

    print_orig(i, j, k, l);
    print_index(i + 1, j - 1 , k + 1, l - 1);
}

// To visualize the data dependecy...
void Foldalign::print_mb_dep(int i, int j, int k, int l, int m, int n) const
{
    return;
    std::cout << i << " " << j << " " << k << " " << l << ":\t"
        << i << " " << m << " " << k << " " << n
        << " + "
        << m + 1 << " " << j << " " << n + 1 << " " << l << "\n";
}

void Foldalign::expand_pos(const int &i, const int &j, const int &k, const int &l)
{
    int score = 0;

    print_score_dep(i, j, k, l);

    score = std::max(score, dp_matrix[index(i + 1, j, k, l)] + Cost::gap);
    score = std::max(score, dp_matrix[index(i, j, k + 1, l)] + Cost::gap);
    score = std::max(score, dp_matrix[index(i, j - 1, k, l)] + Cost::gap);
    score = std::max(score, dp_matrix[index(i, j, k, l - 1)] + Cost::gap);
    score = std::max(score, dp_matrix[index(i + 1, j, k + 1, l)] + Cost::match_score(s1[i], s2[k]));
    score = std::max(score, dp_matrix[index(i, j - 1, k, l - 1)] + Cost::match_score(s1[j], s2[l]));
    score = std::max(score, dp_matrix[index(i + 1, j - 1, k, l)] + Cost::base_score(s1[i], s1[j]) + Cost::gap * 2);
    score = std::max(score, dp_matrix[index(i, j, k + 1, l - 1)] + Cost::base_score(s2[k], s2[l]) + Cost::gap * 2);
    score = std::max(score, dp_matrix[index(i + 1, j - 1, k + 1, l - 1)] +
            Cost::base_score(s1[i], s1[j]) + Cost::base_score(s2[k], s2[l]) +
            Cost::compensation_score(s1[i], s1[j], s2[k], s2[l]));

    for (int m = i + 1; m < j; ++m)
    {
        for (int n = k + 1; n < l; ++n)
        {
            print_mb_dep(i, j, k, l, m, n);
            score = std::max(score, dp_matrix[index(i, m, k, n)] + dp_matrix[index(m + 1, j, n + 1, l)]);
        } //n
    } //m
    if (score > 0)
        dp_matrix[index(i, j, k, l)] = score;
}

void Foldalign::expand_inner_matrix(const int &i, const int &k)
{
    for (int j = i; j < s1_l; ++j)
    {
        for (int l = k; l < s2_l; ++l)
        {
            expand_pos(i, j, k, l);
        } //l
    } //j
}

void Foldalign::expand_inner_matrix_diag(const int &i, const int &k)
{
    for (int inner_diag = i; inner_diag < s1_l; ++inner_diag)
    {
        int j = inner_diag;
        for (int l = k; l < s2_l && j >= i && j < s1_l; ++l, --j)
        {
            expand_pos(i, j, k, l);
        }
    }
    for (int inner_diag = k + 1; inner_diag < s2_l; ++inner_diag)
    {
        int j = s1_l - 1;
        for (int l = inner_diag; l < s1_l && j >= i; --j, ++l)
        {
            expand_pos(i, j, k, l);
        }
    }
}

int Foldalign::sankoff()
{
    std::cout << "Sankoff:"
        << "\nseq1:\t" << s1
        << "\nseq2:\t" << s2
        << "\n";

    for (int i = s1_l - 1; i >= 0; --i)
    {
        for (int k = s2_l - 1; k >= 0; --k)
        {
            expand_inner_matrix(i, k);
        } //k
    } //i
    std::cout << dp_matrix[index(0, s1_l - 1, 0, s2_l - 1)] << std::endl;
    return 0;
}

int Foldalign::diag_sankoff()
{
    std::cout << "Sankoff:"
        << "\nseq1:\t" << s1
        << "\nseq2:\t" << s2
        << "\n";

    for (int outer_diag = 0; outer_diag <= s2_l - 1; ++outer_diag)
    {
        int i = s1_l - 1;
        for (int k = s2_l - 1 - outer_diag; k <= s2_l - 1; ++k, --i)
        {
            expand_inner_matrix_diag(i, k);
        } //i
    } //outer_diag
    for (int outer_diag = s1_l - 2; outer_diag >= 0 ; --outer_diag)
    {
        int i = outer_diag;
        for (int k = 0; k <= s2_l - 1 && i >= 0; ++k, --i)
        {
            expand_inner_matrix_diag(i, k);
        }
    } //outer_diag
    std::cout << dp_matrix[index(0, s1_l - 1, 0, s2_l - 1)] << std::endl;
    return 0;
}

int Foldalign::fold_align()
{
    std::cout << "Foldalign:"
        << "\n\u03B4 = " << delta
        << "\n\u03BB = " << lambda
        << "\nseq1:\t" << s1
        << "\nseq2:\t" << s2
        << "\n";

    for (int i = s1_l - 1; i >= 0; --i)
    {
        for (int k = s2_l - 1; k >= 0; --k)
        {
            int j_end = i + 1 + lambda;
            if (j_end > s1_l)
                j_end = s1_l;

            for (int j = i + 1; j < j_end; ++j)
            {
                int l_begin = j - delta;
                if (l_begin < k + 1)
                    l_begin = k + 1;
                int l_end = j + delta;
                if (l_end > s2_l)
                    l_end = s2_l;

                for (int l = l_begin; l < l_end; ++l)
                {
                    int score = 0;

                    print_score_dep(i, j, k, l);

                    score = std::max(score, dp_matrix[index(i + 1, j, k, l)] + Cost::gap);
                    score = std::max(score, dp_matrix[index(i, j, k + 1, l)] + Cost::gap);
                    score = std::max(score, dp_matrix[index(i, j - 1, k, l)] + Cost::gap);
                    score = std::max(score, dp_matrix[index(i, j, k, l - 1)] + Cost::gap);
                    score = std::max(score, dp_matrix[index(i + 1, j, k + 1, l)] + Cost::match_score(s1[i], s2[k]));
                    score = std::max(score, dp_matrix[index(i, j - 1, k, l - 1)] + Cost::match_score(s1[j], s2[l]));
                    score = std::max(score, dp_matrix[index(i + 1, j - 1, k, l)] + Cost::base_score(s1[i], s1[j]) + Cost::gap * 2);
                    score = std::max(score, dp_matrix[index(i, j, k + 1, l - 1)] + Cost::base_score(s2[k], s2[l]) + Cost::gap * 2);
                    score = std::max(score, dp_matrix[index(i + 1, j - 1, k + 1, l - 1)] +
                                                        Cost::base_score(s1[i], s1[j]) + Cost::base_score(s2[k], s2[l]) +
                                                        Cost::compensation_score(s1[i], s1[j], s2[k], s2[l]));

                    int m_end = i + 1 + lambda;
                    if (m_end > j)
                        m_end = j;
                    for (int m = i + 1; m < m_end; ++m)
                    {
                        int n_begin = m - delta;
                        if (n_begin < k + 1)
                            n_begin = k + 1;
                        int n_end = m + delta;
                        if (n_end > l)
                            n_end = l;

                        for (int n = n_begin; n < n_end; ++n)
                        {
                            print_mb_dep(i, j, k, l, m, n);
                            score = std::max(score, dp_matrix[index(i, m, k, n)] + dp_matrix[index(m + 1, j, n + 1, l)]);
                        } //n
                    } //m
                    if (score > 0)
                        dp_matrix[index(i, j, k, l)] = score;
                } //l
            } //j
        } //k
    } //i
    std::cout << dp_matrix[index(0, s1_l - 1, 0, s2_l - 1)] << std::endl;
    return 0;
}
