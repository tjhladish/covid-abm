#include <unistd.h>
#include "AbcSmc.h"
#include "simulator.h"
#include <cstdlib>
#include "CCRC32.h"
#include "Utility.h"
#include <unordered_set>

using namespace std;

using covid::util::to_string;
using covid::util::mean;
using covid::util::stdev;
using covid::util::max_element;

time_t GLOBAL_START_TIME;

string calculate_process_id(vector<double> &args, string &argstring);
const string SIM_POP = "escambia";
//const string SIM_POP = "florida";
const string HOME_DIR(std::getenv("HOME"));
const string pop_dir = HOME_DIR + "/work/covid-abm/pop-" + SIM_POP;
const string output_dir("/ufrc/longini/tjhladish/");
//const string imm_dir(output_dir + "");

const int RESTART_BURNIN    = 0;
const int FORECAST_DURATION = 500;
const bool RUN_FORECAST     = true;
const int TOTAL_DURATION    = RUN_FORECAST ? RESTART_BURNIN + FORECAST_DURATION : RESTART_BURNIN;
const int JULIAN_TALLY_DATE = 146; // intervention julian date - 1

//Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed) {
Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const string /*process_id*/) {
    Parameters* par = new Parameters();
    par->define_defaults();
    const MmodsScenario ms = (MmodsScenario) args[0];
    assert(ms <= NUM_OF_MMODS_SCENARIOS);
    if (ms < NUM_OF_MMODS_SCENARIOS) { par->numSurveilledPeople = 10e5; } // default is INT_MAX
    par->serial = serial;

    const float T = 0.25;

    par->household_transmissibility = T;
    par->workplace_transmissibility = T;
    par->social_transmissibility    = T;
    //hospitalizedFraction = {0.0, 0.15, 0.9};
    //par->reportedFraction = {0.0, 0.01, 0.5, 0.8, 1.0};      // fraction of asymptomatic, mild, severe, critical, and deaths reported
    par->reportedFraction = {0.0, 0.2, 0.8, 0.8, 1.0};      // fraction of asymptomatic, mild, severe, critical, and deaths reported
    par->randomseed              = rng_seed;
    par->dailyOutput             = false; // turn on for daily prevalence figure, probably uncomment filter in simulator.h for daily output to get only rel. days
    par->periodicOutput          = false;
    par->periodicOutputInterval  = 20;
    par->weeklyOutput            = false;
    par->monthlyOutput           = true;
    par->yearlyOutput            = true;
    par->abcVerbose              = false; // needs to be false to get WHO daily output
    par->runLength               = TOTAL_DURATION;
    par->startDayOfYear          = 1;
    par->annualIntroductionsCoef = 1;

    par->pathogenicityModel = ORIGINAL_LOGISTIC;

    par->traceContacts = false;
    par->mmodsScenario = ms;      // MMODS_CLOSED, MMODS_2WEEKS, MMODS_1PERCENT, MMODS_OPEN

    // These are only initial values for time-structured interventions.  They can be changed dynamically.
//    par->timedInterventions[SCHOOL_CLOSURE].resize(par->runLength, 0.0);
//    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(par->runLength, 0.0);
//    par->timedInterventions[SOCIAL_DISTANCING].resize(par->runLength, 0.0);

    par->timedInterventions[SCHOOL_CLOSURE].resize(30, 0.0);
    par->timedInterventions[SCHOOL_CLOSURE].resize(par->runLength, 1.0);

    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(30, 0.0);
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(par->runLength, 1.0);

    par->timedInterventions[SOCIAL_DISTANCING].resize(30, 0.0);
    par->timedInterventions[SOCIAL_DISTANCING].resize(par->runLength, 0.5);

    par->reportingLag = 14;
    par->symptomToTestLag = 2;

    par->daysImmune = 730;
    par->VES = 0.0;

    //par->hospitalizedFraction = 0.25; // fraction of cases assumed to be hospitalized

    par->probDailyExposure = {1.0e-05};

    par->populationFilename       = pop_dir    + "/population-"         + SIM_POP + ".txt";
    par->locationFilename         = pop_dir    + "/locations-"          + SIM_POP + ".txt";
    par->networkFilename          = pop_dir    + "/network-"            + SIM_POP + ".txt";

    return par;
}

// Take a list of values, return original indices sorted by value
vector<int> ordered(vector<int> const& values) {

    vector<pair<int,int> > pairs(values.size());
    for(unsigned int pos=0; pos<values.size(); pos++) {
        pairs[pos] = make_pair(values[pos],pos);
    }

    //bool comparator ( const mypair& l, const mypair& r) { return l.first < r.first; }
    std::sort( pairs.rbegin(), pairs.rend() ); // sort greatest to least
    vector<int> indices(values.size());
    for(unsigned int i=0; i < pairs.size(); i++) indices[i] = pairs[i].second;

    return indices;
}


string calculate_process_id(vector<double> &args, string &argstring) {
    // CCRC32 checksum based on string version of argument values
    CCRC32 crc32;
    crc32.Initialize();

    for (unsigned int i = 0; i < args.size(); i++) argstring += to_string((double) args[i]) + " ";

    const unsigned char* argchars = reinterpret_cast<const unsigned char*> (argstring.c_str());
    const int len = argstring.length();
    const int process_id = crc32.FullCRC(argchars, len);

    return to_string(process_id);
}


string report_process_id (vector<double> &args, const unsigned long int serial, const time_t start_time) {
    double dif = difftime (start_time, GLOBAL_START_TIME);

    string argstring;
    const string process_id = calculate_process_id(args, argstring);

    cerr << "pid in report_process_id (num args = " << args.size() << "): " << process_id << endl;
    stringstream ss;
    ss << "begin " << process_id << " " << dec << serial << " " << dif << " " << argstring << endl;
    string output = ss.str();
    fputs(output.c_str(), stderr);

    return to_string(process_id);
}


void append_if_finite(vector<double> &vec, double val) {
    if (isfinite(val)) {
        vec.push_back((double) val);
    } else {
        vec.push_back(0);
    }
}

int julian_to_sim_day (const Parameters* par, const int julian, const int intervention_year) {
    int startDate = intervention_year*365 + julian - par->startDayOfYear;
    if (julian < par->startDayOfYear) { // start intervention in following year
        startDate += 365;
    }
    return startDate;
}


vector<double> tally_counts(const Parameters* par, Community* community, int pre_intervention_output) {
    const int discard_days = julian_to_sim_day(par, JULIAN_TALLY_DATE, RESTART_BURNIN-pre_intervention_output);

    //vector< vector<int> > severe      = community->getNumSevereCases();
    vector<size_t> symptomatic = community->getNumNewlySymptomatic();

    //vector< vector<int> > infected    = community->getNumNewlyInfected();
    const int num_years = FORECAST_DURATION + pre_intervention_output - 1; // typically 55
    vector<size_t> s_tally(num_years+1, 0);

    vector<double> metrics(num_years, 0.0);
    for (int t=discard_days; t<par->runLength; t++) {
        // use epidemic years, instead of calendar years
        const int y = (t-discard_days)/365;
        s_tally[y] += symptomatic[t];
        //cout << "d,i:" << t << "," << infected[0][t] + infected[1][t] + infected[2][t] + infected[3][t] << endl;
    }

    for (int y = 0; y<num_years; ++y) {
        metrics[y] += s_tally[y];
    }

    return metrics;
}

vector<double> simulator(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par* mp = nullptr) {
    gsl_rng_set(RNG, rng_seed); // seed the rng using sys time and the process id
    // initialize bookkeeping for run
    time_t start, end;
    time (&start);
    //const string process_id = report_process_id(args, serial, mp, start);
    //vector<double> abc_args(&args[0], &args[8]);
    vector<double> abc_args(args);
    //const unsigned int realization = 0; //(int) args[9];

    //const string process_id = report_process_id(abc_args, serial, start) + "." + to_string(realization);
    const string process_id = to_string(rng_seed);
    report_process_id(args, serial, start);

    cerr << "SCENARIO " << rng_seed;
    for (auto _p: args) { cerr << " " << _p; } cerr << endl;

    Parameters* par = define_simulator_parameters(args, rng_seed, serial, process_id);
    //00   "mild_rf",
    //02   "severe_rf",
    //03   "sec_path",
    //04   "sec_sev",
    //05   "pss_ratio",
    //06   "exp_coef",
    //08   "realization",

    Community* community = build_community(par);

    seed_epidemic(par, community);
    simulate_epidemic(par, community, process_id);
    //community->getMeanNumSecondaryInfections();
    time (&end);
    double dif = difftime (end,start);

//    const int pre_intervention_output = 5; // years
//    const int desired_intervention_output = FORECAST_DURATION - 1;
    vector<double> metrics(0);// = tally_counts(par, community, pre_intervention_output);

    const vector<size_t> infections       = community->getNumNewlyInfected();
    const vector<size_t> reported_cases   = community->getNumDetectedCasesReport();
    const vector<size_t> hospitalizations = community->getNumHospPrev();
//    const vector<size_t> icu              = community->getNumIcuInc();
    const vector<size_t> deaths           = community->getNumDetectedDeaths();
    const vector<float> lockdown          = community->getTimedIntervention(NONESSENTIAL_BUSINESS_CLOSURE);

//    string header = "scenario replicate day infections deaths hosp_prev detected ne_closed";
//    const int scenario = 0;
//    const size_t replicate = rng_seed;
//    for (int day = 0; day < par->runLength; ++day) {
//        cerr
//            << scenario
//            << " " << replicate
//            << " " << day
//            << " " << infections[day]
//            << " " << deaths[day]
//            << " " << hospitalizations[day]
//            << " " << reported_cases[day]
//            << " " << lockdown[day]
//            << endl;
//    }

    stringstream ss;
    ss << mp->mpi_rank << " end " << hex << process_id << " " << dec << dif << " ";

    for (auto i: args) ss << i << " ";
    for (auto i: metrics) ss << i << " ";
    ss << endl;

    string output = ss.str();
    fputs(output.c_str(), stderr);

    delete par;
    delete community;
    return metrics;
}

/*void usage() {
    cerr << "\n\tUsage: ./abc_sql abc_config_sql.json --process\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate -n <number of simulations per database write>\n\n";

}*/


int main(int argc, char* argv[]) {
    assert(argc == 3);
    vector<double> sim_args = {atof(argv[1])};
    size_t seed = atoi(argv[2]);
    size_t serial = 0;

    simulator(sim_args, seed, serial);

    return 0;
}
