#ifndef _SANKOFF_GPU_H
#define _SANKOFF_GPU_H
#include <string>

long long int dp_matrix_calc_total_size(long long int s1, long long int s2);
int dp_matrix_calc_delta(int i, int j, int k, int l, const int &s1_l, const int &s2_l);
bool dp_matrix_check_border(const int &i, const int &j, const int &k, const int &l, const int &s1_l, const int &s2_l);
int dp_matrix_get_pos(int *dp_matrix, const int &i, const int &j, const int &k, const int &l, const int &s1_l, const int &s2_l);
void dp_matrix_put_pos(int *dp_matrix, const int &i, const int &j, const int &k, const int &l, const int &val, const int &s1_l, const int &s2_l);

class Sankoff_GPU {
    public:
        Sankoff_GPU(const std::string &seq1, const std::string &seq2);
        virtual ~Sankoff_GPU();
        int diag_sankoff(); //Run a pure sankoff algorithm

    private:
        void expand_inner_matrix_diag(const int &i, const int &k);
        void expand_pos(int *dp_matrix, const int &i, const int &j, const int &k, const int &l, const int &s1_l, const int &s2_l);

        //device_members
        int *dp_matrix;
        std::string s1; 
        int s1_l;
        std::string s2;
        int s2_l;
};
#endif
