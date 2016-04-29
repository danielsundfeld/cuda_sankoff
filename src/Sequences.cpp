/*!
 * \class Sequences
 * \author Daniel Sundfeld
 * \copyright MIT License
 */
#include "Sequences.h"

#include <iostream>
#include <string>

Sequences::Sequences()
{
}

//! Number of sequences
int Sequences::n_seq = 0;
//! Singleton instance
Sequences Sequences::instance;
//! Destination coord

//! Save the string \a x as an Sequence
int Sequences::set_seq(const std::string &x)
{
    seqs.push_back(x);
    ++n_seq;
    return n_seq;
}
