#include "Foldalign.h"

#include <iostream>

#include "Cost.h"

Foldalign::Foldalign(const std::string &seq1, const std::string &seq2,
                     const int lmbd, const int dlt)
: Sankoff(seq1, seq2)
{
    lambda = lmbd;
    delta = dlt;
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
                    float score = 0;

                    print_score_dep(i, j, k, l);

                    score = std::max(score, dp_matrix.get_pos(i + 1, j, k, l) + Cost::gap);
                    score = std::max(score, dp_matrix.get_pos(i, j, k + 1, l) + Cost::gap);
                    score = std::max(score, dp_matrix.get_pos(i, j - 1, k, l) + Cost::gap);
                    score = std::max(score, dp_matrix.get_pos(i, j, k, l - 1) + Cost::gap);
                    score = std::max(score, dp_matrix.get_pos(i + 1, j, k + 1, l) + Cost::match_score(s1[i], s2[k]));
                    score = std::max(score, dp_matrix.get_pos(i, j - 1, k, l - 1) + Cost::match_score(s1[j], s2[l]));
                    score = std::max(score, dp_matrix.get_pos(i + 1, j - 1, k, l) + Cost::base_score(s1[i], s1[j]) + Cost::gap * 2);
                    score = std::max(score, dp_matrix.get_pos(i, j, k + 1, l - 1) + Cost::base_score(s2[k], s2[l]) + Cost::gap * 2);
                    score = std::max(score, dp_matrix.get_pos(i + 1, j - 1, k + 1, l - 1) +
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
                            score = std::max(score, dp_matrix.get_pos(i, m, k, n) + dp_matrix.get_pos(m + 1, j, n + 1, l));
                        } //n
                    } //m
                    if (score > 0)
                        dp_matrix.put_pos(i, j, k, l, score);
                } //l
            } //j
        } //k
    } //i
    std::cout << dp_matrix.get_pos(0, s1_l - 1, 0, s2_l - 1) << std::endl;
    return 0;
}
