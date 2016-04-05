#include "Sankoff_GPU.h"

#include <iostream>

#include "Cost.h"

long long int dp_matrix_calc_total_size(long long int s1, long long int s2)
{
    return (((1 + s1) * s1) / 2 ) * (((1 + s2) * s2) / 2);
}

int dp_matrix_calc_delta(int i, int j, int k, int l, const int &s1_l, const int &s2_l)
{
    i = s1_l - i;
    j = s1_l - j;
    k = s2_l - k;
    l = s2_l - l;

    int delta_i = ((1 + s1_l) * s1_l / 2) * (i * (i - 1) / 2);
    int delta_k = i * (k * (k - 1) / 2);
    int delta_mi = (j - 1) * k + l - 1;
    return delta_i + delta_k + delta_mi;
}

bool dp_matrix_check_border(const int &i, const int &j, const int &k, const int &l, const int &s1_l, const int &s2_l)
{
    if (j < i)
        return false;
    if (l < k)
        return false;
    if (i < 0 || j < 0 || k < 0 || l < 0)
        return false;
    if (i > s1_l || j > s1_l || k > s2_l || l > s2_l)
        return false;
    return true;
}

int dp_matrix_get_pos(int *dp_matrix, const int &i, const int &j, const int &k, const int &l, const int &s1_l, const int &s2_l)
{
    if (dp_matrix_check_border(i, j, k, l, s1_l, s2_l) == false)
        return -1024;

    return dp_matrix[dp_matrix_calc_delta(i, j, k, l, s1_l, s2_l)];
}

void dp_matrix_put_pos(int *dp_matrix, const int &i, const int &j, const int &k, const int &l, const int &val, const int &s1_l, const int &s2_l)
{
    dp_matrix[dp_matrix_calc_delta(i, j, k, l, s1_l, s2_l)] = val;
}

Sankoff_GPU::Sankoff_GPU(const std::string &seq1, const std::string &seq2)
{
    s1 = seq1;
    s1_l = (int) s1.length(); //TODO: throw exception if negative
    s2 = seq2;
    s2_l = (int) s2.length();

    cudaMalloc(&dp_matrix, dp_matrix_calc_total_size(s1_l, s2_l) * sizeof(int));
}

Sankoff_GPU::~Sankoff_GPU()
{
    cudaFree(dp_matrix);
}

void Sankoff_GPU::expand_pos(int *dp_matrix, const int &i, const int &j, const int &k, const int &l, const int &s1_l, const int &s2_l)
{
    int score = 0;

    score = std::max(score, dp_matrix_get_pos(dp_matrix, i + 1, j, k, l, s1_l, s2_l) + Cost::gap);
    score = std::max(score, dp_matrix_get_pos(dp_matrix, i, j, k + 1, l, s1_l, s2_l) + Cost::gap);
    score = std::max(score, dp_matrix_get_pos(dp_matrix, i, j - 1, k, l, s1_l, s2_l) + Cost::gap);
    score = std::max(score, dp_matrix_get_pos(dp_matrix, i, j, k, l - 1, s1_l, s2_l) + Cost::gap);
    score = std::max(score, dp_matrix_get_pos(dp_matrix, i + 1, j, k + 1, l, s1_l, s2_l) + Cost::match_score(s1[i], s2[k]));
    score = std::max(score, dp_matrix_get_pos(dp_matrix, i, j - 1, k, l - 1, s1_l, s2_l) + Cost::match_score(s1[j], s2[l]));
    score = std::max(score, dp_matrix_get_pos(dp_matrix, i + 1, j - 1, k, l, s1_l, s2_l) + Cost::base_score(s1[i], s1[j]) + Cost::gap * 2);
    score = std::max(score, dp_matrix_get_pos(dp_matrix, i, j, k + 1, l - 1, s1_l, s2_l) + Cost::base_score(s2[k], s2[l]) + Cost::gap * 2);
    score = std::max(score, dp_matrix_get_pos(dp_matrix, i + 1, j - 1, k + 1, l - 1, s1_l, s2_l) +
            Cost::base_score(s1[i], s1[j]) + Cost::base_score(s2[k], s2[l]) +
            Cost::compensation_score(s1[i], s1[j], s2[k], s2[l]));

    for (int m = i + 1; m < j; ++m)
    {
        for (int n = k + 1; n < l; ++n)
        {
            score = std::max(score, dp_matrix_get_pos(dp_matrix, i, m, k, n, s1_l, s2_l) + dp_matrix_get_pos(dp_matrix, m + 1, j, n + 1, l, s1_l, s2_l));
        } //n
    } //m
    dp_matrix_put_pos(dp_matrix, i, j, k, l, score, s1_l, s2_l);
}

void Sankoff_GPU::expand_inner_matrix_diag(const int &i, const int &k)
{
    for (int inner_diag = i; inner_diag < s1_l; ++inner_diag)
    {
#pragma omp parallel for schedule(dynamic,1)
        for (int l = k; l < s2_l; ++l)
        {
            int j = inner_diag - (l - k);
            if (j >= i && j < s1_l)
                expand_pos(dp_matrix, i, j, k, l, s1_l, s2_l);
        }
    }
    for (int inner_diag = k + 1; inner_diag < s2_l; ++inner_diag)
    {
#pragma omp parallel for schedule(dynamic,1)
        for (int l = inner_diag; l < s1_l; ++l)
        {
            int j = s1_l - 1 - (l - inner_diag);
            if (j >= i)
                expand_pos(dp_matrix, i, j, k, l, s1_l, s2_l);
        }
    }
}

int Sankoff_GPU::diag_sankoff()
{
    std::cout << "Sankoff_GPU:"
        << "\nseq1:\t" << s1
        << "\nseq2:\t" << s2
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
    //TODO copy to CPU
    //std::cout << dp_matrix_get_pos(0, s1_l - 1, 0, s2_l - 1) << std::endl;
    std::cout << "FIM\n";
    return 0;
}
