#include <chrono>
#include <unistd.h>
#include "AbcSmc.h"
#include <cstdlib>
#include "CCRC32.h"
#include "Utility.h"
#include <math.h>

#if __has_include("local.h")
#include "local.h"
#endif

using namespace std;

using covid::util::to_string;
using covid::util::mean;
using covid::util::stdev;
using covid::util::max_element;

time_t GLOBAL_START_TIME;

string calculate_process_id(vector<double> &args, string &argstring) {
    // CCRC32 checksum based on string version of argument values
    CCRC32 crc32;
    crc32.Initialize();

    for (size_t i = 0; i < args.size(); i++) argstring += to_string((double) args[i]) + " ";

    const unsigned char* argchars = reinterpret_cast<const unsigned char*> (argstring.c_str());
    const int len = argstring.length();
    const int process_id = crc32.FullCRC(argchars, len);

    return to_string(process_id);
}

#ifdef LOCAL_HEADER
const std::string HOME_DIR = getHOME();
#else
const std::string HOME_DIR = std::getenv("HOME");
#endif

int main(int argc, char* argv[]) {
    if (not (argc == 1) ) {
        cerr << "ERROR: please provide JSON config file" << endl;
        exit(100);
    }

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(string(argv[1]));

    return 0;
}
