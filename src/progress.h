#ifndef PROGRESS
#define PROGRESS

#include <time.h>
#include <string>

class Progress {
private:
    int end;
    std::string name;
    int freq;
    int prog;
    int last;
    time_t lasttime;
    int maxwait;
    
    void createBar ();
public:
    Progress (int e, std::string n);
    Progress (int e, std::string n, int f);
    void increment ();
    void start ();
    void finish ();
};

#endif