#ifndef _DPMATRIX_H
#define _DPMATRIX_H
class DPMatrix {
    public:
        DPMatrix(const int &s1_l, const int &s2_l);
        ~DPMatrix();
        int get_pos(const int &i, const int &j, const int &k, const int &l);
        void put_pos(const int &i, const int &j, const int &k, const int &l, const int &val);
        long long int get_total_size() const;

    private:
        int calc_delta(int i, int j, int k, int l) const;
        long long int calc_total_size(long long int s1, long long int s2) const;
        bool check_border(const int &i, const int &j, const int &k, const int &l) const;

        int *dp_matrix;
        const int s1_l;
        const int s2_l;
};
#endif
