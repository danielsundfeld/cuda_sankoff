/*!
 * \class TimeCounter
 * \author Daniel Sundfeld
 * \copyright MIT License
 *
 * \brief Class that count elapsed time between the creation and destruction
 */
#ifndef _TIMECOUNTER_H
#define _TIMECOUNTER_H
#include <sys/time.h>
#include <iostream>
#include <string>

class TimeCounter
{
    public:
        TimeCounter(const std::string &msg = "");
        ~TimeCounter();

    private:
        struct timeval m_begin;
        std::string m_msg;
};
#endif
