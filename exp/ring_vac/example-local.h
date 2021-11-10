#include <string>

#ifndef LOCAL_HEADER
#define LOCAL_HEADER
//CHANGE LOCAL PATHS BELOW
std::string localHome = "LOCAL PATH GOES HERE";

// DO NOT ALTER FUNCTION BELOW
const std::string getHOME() {
    if(localHome.size() != 0) {
        return localHome;
    }
    else {
        return std::getenv("HOME");
    }
}
#endif
