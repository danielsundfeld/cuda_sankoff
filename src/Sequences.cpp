/*!
 * \class Sequences
 * \author Daniel Sundfeld
 * \copyright MIT License
 */
#include "Sequences.h"

#include <iostream>
#include <string>

//! Number of sequences
int Sequences::n_seq = 0;
//! Singleton instance
Sequences Sequences::instance;

//! Save the string \a x as an Sequence
int Sequences::set_seq(const std::string &x)
{
    if (n_seq && seqs[0].length() < x.length())
        seqs.insert(seqs.begin(), x);
    else
        seqs.push_back(x);
    ++n_seq;
    return n_seq;
}
