#ifndef _DPMATRIX_GPU_H
#define _DPMATRIX_GPU_H

#define MAX_SEQ_SIZE 1024

struct sequences {
    char s1[MAX_SEQ_SIZE];
    int s1_l;
    char s2[MAX_SEQ_SIZE];
    int s2_l;
};

long long int dp_matrix_calc_total_size(long long int s1, long long int s2);
int dp_matrix_calc_delta(int i, int j, int k, int l, const int &s1_l, const int &s2_l);
bool dp_matrix_check_border(const int &i, const int &j, const int &k, const int &l, const int &s1_l, const int &s2_l);
int dp_matrix_get_pos(int *dp_matrix, const int &i, const int &j, const int &k, const int &l, const int &s1_l, const int &s2_l);
void dp_matrix_put_pos(int *dp_matrix, const int &i, const int &j, const int &k, const int &l, const int &val, const int &s1_l, const int &s2_l);
#endif
