/*
 *  cp_time.h
 *  cross platform time
 *
 *  Created by Jimmy on 12/14/10.
 *  Copyright 2010 OSU. All rights reserved.
 *
 */
#ifndef CP_TIMER_H
#define CP_TIMER_H
#include <string>
#include <stdio.h>

#ifdef WIN32
#include <Windows.h>
#include <mmsystem.h> //Winmm.lib
//typedef  long long __int64;
#else
#include <sys/time.h>
#endif

class Timer
{
    long long time1, time2; // in  us
    static long long getTime() {
#ifdef WIN32
        return (__int64)timeGetTime() * 1000;
#else
        struct timezone Tzp;
        struct timeval Tp;
        int stat = gettimeofday (&Tp, &Tzp);
        if (stat) printf("Error return from gettimeofday: %d\n",stat);
        return ((long long)Tp.tv_sec*1000000 + Tp.tv_usec);        
#endif
    }
    
public:

    Timer() {
    }
    
    inline void start() {
        time1 = getTime();
    }
    
    inline void end() {
        time2 = getTime();
    }
    
    static inline long long getTimeMS()
    {
        return getTime()/1000;
    }
    
    static inline long long getTimeUS()
    {
        return getTime();
    }
    
    static inline struct tm *getTimeInfo()
    {
        time_t rawtime;
        struct tm * info;
        time ( &rawtime );
        info = localtime ( &rawtime );
        return info;
    }
    static inline std::string getTimeString()
    {
        struct tm *info = getTimeInfo();
        char s[256];
        sprintf ( s, "%02d%02d%02d_%02d%02d%02d", info->tm_year-100, info->tm_mon, info->tm_mday, info->tm_hour, info->tm_min, info->tm_sec );
        return std::string(s);
    }

    
    
        
    inline long long getElapsedMS()
    {
        return (time2-time1)/1000;
    }
    
    inline long long getElapsedUS()
    {
        return time2-time1;
    }
    inline long long getStartTimeMS()
    {
        return time2/1000;
    }
    inline long long getEndTimeMS()
    {
        return time1/1000;
    }};
#endif 