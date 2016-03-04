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
    bool show;
    
    void createBar ();
public:
    Progress (int e, std::string n, bool s);
    Progress (int e, std::string n, int f, bool s);
    void increment ();
    void start ();
    void finish ();
};

#endif