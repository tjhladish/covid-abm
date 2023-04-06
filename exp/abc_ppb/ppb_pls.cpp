#include <chrono>
#include <unistd.h>
#include "AbcSmc.h"
#include <cstdlib>
#include "CCRC32.h"
#include <math.h>

#if __has_include("local.h")
#include "local.h"
#endif

#ifdef LOCAL_HEADER
const std::string HOME_DIR = getHOME();
#else
const std::string HOME_DIR = std::getenv("HOME");
#endif

using namespace std;

int main(int argc, char* argv[]) {
    if (not (argc == 2) ) {
        cerr << "ERROR: please provide JSON config file" << endl;
        exit(100);
    }

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(string(argv[1]));

    sqdb::Db db(abc->get_database_filename().c_str());
    vector< vector<int> > serials;
    abc->read_SMC_sets_from_database (db, serials);
    //_filter_particles( t, _particle_metrics[t], _particle_parameters[t], next_pred_prior_size );
    const int set = 0;             // dummy value in this case
    ABC::Mat2D parameters = abc->get_particle_parameters()[set];
    ABC::Mat2D metrics    = abc->get_particle_metrics()[set];

    const int pred_prior_size = 0; // dummy value in this case
    PLS_Model plsm = abc->run_PLS(parameters, metrics, parameters.rows(), parameters.cols());
    //PLS_Model plsm = abc->_filter_particles(set,  particle_parameters[set], particle_metrics[set],  pred_prior_size);

    return 0;
}
