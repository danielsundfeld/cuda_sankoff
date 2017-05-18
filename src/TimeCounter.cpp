/*!
 * \class TimeCounter
 * \author Daniel Sundfeld
 * \copyright GPL License
 */
#include "TimeCounter.h"

#include <string>
#include <stdio.h>
#include <unistd.h>

TimeCounter::TimeCounter(const std::string &msg)
: m_msg(msg)
{
    gettimeofday(&m_begin, NULL);
}

// 
int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y)
{
    /* Perform the carry for the later subtraction by updating y. */
    if (x->tv_usec < y->tv_usec) {
        int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
        y->tv_usec -= 1000000 * nsec;
        y->tv_sec += nsec;
    }
    if (x->tv_usec - y->tv_usec > 1000000) {
        int nsec = (x->tv_usec - y->tv_usec) / 1000000;
        y->tv_usec += 1000000 * nsec;
        y->tv_sec -= nsec;
    }

    /* Compute the time remaining to wait.
       tv_usec is certainly positive. */
    result->tv_sec = x->tv_sec - y->tv_sec;
    result->tv_usec = x->tv_usec - y->tv_usec;

    /* Return 1 if result is negative. */
    return x->tv_sec < y->tv_sec;
}

TimeCounter::~TimeCounter()
{
    struct timeval end, result;

    gettimeofday(&end, NULL);
    timeval_subtract(&result, &end, &m_begin);
    printf("%s: %ld:%02ld.%06ld\n", m_msg.c_str(), result.tv_sec / 60, result.tv_sec % 60, result.tv_usec);
}
