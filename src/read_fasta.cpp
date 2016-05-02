/*!
 * \author Daniel Sundfeld
 * \copyright MIT License
 */
#include "read_fasta.h"

#include <fstream>
#include <iostream>
#include <string>

#include "Sequences.h"

int read_fasta_file_core(const std::string &name)
{
    std::ifstream file(name.c_str());
    Sequences *sequences = Sequences::get_instance();

    if (!file.is_open())
    {
        std::cout << "Can't open file " << name << std::endl;
        return -1;
    }

    while (!file.eof())
    {
        std::string seq;
        while (!file.eof())
        {
            std::string buf;
            getline(file, buf);
            if ((buf.length() <= 0) || (buf.at(0) == '>'))
                break;
            seq.append(buf);
        }
        if (!seq.empty())
            sequences->set_seq(seq);
    }
    return 0;
}

/*!
 * Read the \a name fasta file, loading it to the Sequences singleton.
 * Fail if the number of sequences is non-null and not equal to \a limit
 */
int read_fasta_file(const std::string &name, const int &check)
{
    try
    {
        int ret = read_fasta_file_core(name);
        if (ret == 0 && check != 0 && Sequences::get_nseq() != check)
        {
            std::cerr << "Invalid fasta file: must have " << check << " sequences.\n";
            return -1;
        }
        return ret;
    }
    catch (std::exception &e)
    {
        std::cerr << "Reading file fatal error: " << e.what() << std::endl;
    }
    catch (...)
    {
        std::cerr << "Unknown fatal error while reading file!\n";
    }
    return -1;
}
