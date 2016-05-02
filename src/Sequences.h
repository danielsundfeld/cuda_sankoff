/*!
 * \class Sequences
 * \author Daniel Sundfeld
 * \copyright MIT License
 *
 * \brief Singleton that holds all sequences being aligned
 */
#ifndef _SEQUENCES_H
#define _SEQUENCES_H
#include <string>
#include <vector>

class Sequences
{
    public:
        static Sequences* get_instance() { return &instance; };
        int set_seq(const std::string &x);
        const std::string& get_seq(int x) const { return seqs.at(x); };
        static int get_nseq() { return n_seq; };

    private:
        static int n_seq;
        static Sequences instance;
        std::vector<std::string> seqs;
        Sequences() { };
};
#endif
