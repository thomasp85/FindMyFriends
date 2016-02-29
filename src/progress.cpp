#include <Rcpp.h>
#include "progress.h"

using namespace Rcpp;

Progress::Progress (int e, std::string n) {
    end = e;
    name = n;
    freq = e / 50;
    prog = 0;
    last = 0;
    std::time(&lasttime);
    maxwait = 10;
}
Progress::Progress (int e, std::string n, int f) {
    end = e;
    name = n;
    if (e/f < 50) {
        freq = e / 50;
    } else {
        freq = f;
    }
    prog = 0;
    last = 0;
    std::time(&lasttime);
    maxwait = 10;
}
void Progress::createBar() {
    int p = 50 * double(prog) / end;
    if (p > 50) p = 50;
    int countWhite = std::to_string(end).size() - std::to_string(prog).size();
    Rcout << "\r" << name << " |" << std::string(p, '=') << std::string(50-p, ' ') << "| " << std::string(countWhite, ' ') << prog << "/" << end << std::flush;
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
    std::time_t timer;
    std::time(&timer);
    
    if (end == prog) {
        createBar();
    } else if (last > freq || maxwait < std::difftime(timer, lasttime)) {
        R_CheckUserInterrupt();
        lasttime = timer;
        last = 0;
        createBar();
    }
}