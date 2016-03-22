#include <time.h>
#include <Rcpp.h>
#include "progress.h"

using namespace Rcpp;

int nDigits(int n) {
    int number_of_digits = 0;
    do {
        ++number_of_digits; 
        n /= 10;
    } while (n);
    return number_of_digits;
}

Progress::Progress (int e, std::string n, bool s) {
    show = s;
    end = e;
    name = n;
    freq = e / 50;
    prog = 0;
    last = 0;
    time(&lasttime);
    maxwait = 10;
}
Progress::Progress (int e, std::string n, int f, bool s) {
    show = s;
    end = e;
    name = n;
    if (e/f < 50) {
        freq = e / 50;
    } else {
        freq = f;
    }
    prog = 0;
    last = 0;
    time(&lasttime);
    maxwait = 10;
}
void Progress::createBar() {
    if (show) {
        int p = 50 * double(prog) / end;
        if (p > 50) p = 50;
        int countWhite = nDigits(end) - nDigits(prog);
        Rcout << "\r" << name << " |" << std::string(p, '=') << std::string(50-p, ' ') << "| " << std::string(countWhite, ' ') << prog << "/" << end << std::flush;
    }
}
void Progress::start() {
    createBar();
}
void Progress::finish() {
    prog = end;
    createBar();
}
void Progress::increment() {
    ++prog;
    if (prog > end) prog = end;
    ++last;
    time_t timer;
    time(&timer);
    
    if (end == prog) {
        createBar();
    } else if (last >= freq || maxwait < difftime(timer, lasttime)) {
        R_CheckUserInterrupt();
        lasttime = timer;
        last = 0;
        createBar();
    }
}