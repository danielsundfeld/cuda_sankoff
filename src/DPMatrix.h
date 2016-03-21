#ifndef _DPMATRIX_H
#define _DPMATRIX_H
#include <map>
#include <string>
#include <tuple>

typedef std::tuple<int, int, int, int> index;

class DPMatrix {
    public:
        DPMatrix(const int &s1_l, const int &s2_l);
        int get_pos(const int &i, const int &j, const int &k, const int &l);
        void put_pos(const int &i, const int &j, const int &k, const int &l, const int &val);

    private:
        int calc_delta(int i, int j, int k, int l) const;
        std::map<index, int> matrix;

        const int s1_l;
        const int s2_l;
};
#endif
