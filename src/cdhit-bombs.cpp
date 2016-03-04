#include <Rcpp.h>
#include "cdhit-common.h"

using namespace Rcpp;

Function rWarning("warning");

void bomb_error(const char *message)
{
    clear_temps();
    std::string error = "\nFatal Error:\n%s\nProgram halted !!\n\n";
    stop(error + message);
} // END void bomb_error

void bomb_error(const char *message, const char *message2)
{
    clear_temps();
    std::string error = "\nFatal Error:\n%s\nProgram halted !!\n\n";
    stop(error + message + message2);
} // END void bomb_error

void bomb_warning(const char *message)
{
    std::string warn = "\nWarning:\n%s\nNot fatal, but may affect results !!\n\n";
    rWarning(warn + message);
} // END void bomb_warning


void bomb_warning(const char *message, const char *message2)
{
    std::string warn = "\nWarning:\n%s\nNot fatal, but may affect results !!\n\n";
    rWarning(warn + message + message2);
} // END void bomb_warning
