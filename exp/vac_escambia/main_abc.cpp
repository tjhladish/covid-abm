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

const int RESTART_BURNIN       = 0;
const int FORECAST_DURATION    = 364;
const bool RUN_FORECAST        = true;
const int TOTAL_DURATION       = RUN_FORECAST ? RESTART_BURNIN + FORECAST_DURATION : RESTART_BURNIN;
const size_t JULIAN_TALLY_DATE = 146; // intervention julian date - 1
const size_t JULIAN_START_YEAR = 2020;

//Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed) {
Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const string /*process_id*/) {
    Parameters* par = new Parameters();
    par->define_defaults();
    par->serial = serial;

    const float T = 0.225;

    par->household_transmissibility = T;
    par->workplace_transmissibility = T;
    par->school_transmissibility    = T;
    par->social_transmissibility    = T/10.0;
    //hospitalizedFraction = {0.0, 0.15, 0.9};
    //par->reportedFraction = {0.0, 0.01, 0.5, 0.8, 1.0};      // fraction of asymptomatic, mild, severe, critical, and deaths reported
    par->reportedFraction = {0.0, 0.2, 0.75, 0.75, 0.75};      // fraction of asymptomatic, mild, severe, critical, and deaths reported
    par->randomseed              = rng_seed;
    par->dailyOutput             = false; // turn on for daily prevalence figure, probably uncomment filter in simulator.h for daily output to get only rel. days
    par->periodicOutput          = false;
    par->periodicOutputInterval  = 20;
    par->weeklyOutput            = false;
    par->monthlyOutput           = true;
    par->yearlyOutput            = true;
    par->abcVerbose              = false; // needs to be false to get WHO daily output
    par->runLength               = TOTAL_DURATION;
    par->julianYear              = JULIAN_START_YEAR;
    par->startDayOfYear          = 45;
//    par->startDayOfYear          = 1;
    par->annualIntroductionsCoef = 1;

    //par->pathogenicityModel = ORIGINAL_LOGISTIC;
    par->traceContacts = true;
    //par->mmodsScenario = ms;      // MMODS_CLOSED, MMODS_2WEEKS, MMODS_1PERCENT, MMODS_OPEN

    // These are only initial values for time-structured interventions.  They can be changed dynamically.
    par->timedInterventions[SCHOOL_CLOSURE].resize(30, 0.0);
    const size_t aug31 = 243 - par->startDayOfYear;
    const size_t school_closed_duration = aug31 < par->runLength ? aug31 : par->runLength;
    par->timedInterventions[SCHOOL_CLOSURE].resize(school_closed_duration, 1.0);
    par->timedInterventions[SCHOOL_CLOSURE].resize(par->runLength, 0.0);

    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].clear();
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(45, 0.0);
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(75, 1.0);
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(par->runLength, 0.0);

    const float mobility_logit_shift   = -2.7;
    const float mobility_logit_stretch = 2;
    // using createSocialDistancingModel() defines timedInterventions[SOCIAL_DISTANCING] values
    // NB: startDayOfYear, runLength, and julianYear must be defined in par before a social distancing model can be meaningfully created!
    // TODO - make that dependency something the user doesn't have to know or think about
    //par->createSocialDistancingModel(pop_dir + "/safegraph_mobility_index.csv", mobility_logit_shift, mobility_logit_stretch);
    // plot(as.Date(d$date), shiftstretch(1-d$smooth_ma7, shift=-2.7, stretch=2), ylim=c(0,1), type='l')
    par->createSocialDistancingModel(pop_dir + "/sgmi_fl_comp.csv", mobility_logit_shift, mobility_logit_stretch);

    //par->defaultReportingLag = 14;
    par->createReportingLagModel(pop_dir + "/case_report_delay.csv");
    par->symptomToTestLag = 2;
    par->deathReportingLag = 4;

    par->daysImmune = 730;
    par->VES = 0.0;

    //par->hospitalizedFraction = 0.25; // fraction of cases assumed to be hospitalized
    par->numInitialInfected = 1;
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


int julian_to_sim_day (const Parameters* par, const size_t julian, const int intervention_year) {
    int startDate = intervention_year*365 + julian - par->startDayOfYear;
    if (julian < par->startDayOfYear) { // start intervention in following year
        startDate += 365;
    }
    return startDate;
}


vector<double> tally_counts(const Parameters* par, Community* community, int discard_days) {
    const size_t num_weeks = (par->runLength - discard_days)/7;

    //vector<size_t> infected    = community->getNumNewlyInfected();
    //vector< vector<int> > severe      = community->getNumSevereCases();
    vector<size_t> symptomatic = community->getNumNewlySymptomatic();

    // pair of num of primary infections starting this day, and mean num secondary infections they cause
    vector<pair<size_t, double>> R = community->getMeanNumSecondaryInfections();
    vector<size_t> Rt_incidence_tally(num_weeks, 0);

    vector<double> metrics(num_weeks*2, 0.0);
    for (size_t t = discard_days; t < par->runLength; t++) {
        const size_t w = (t-discard_days)/7; // which reporting week are we in?
        metrics[w] += symptomatic[t];
        metrics[num_weeks + w] += R[t].first > 0 ? R[t].first*R[t].second : 0;
        Rt_incidence_tally[w] += R[t].first;
    }

    for (size_t w = 0; w < num_weeks; ++w) {
        metrics[num_weeks + w] /= Rt_incidence_tally[w] > 0 ? Rt_incidence_tally[w] : 1.0;
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
    const bool vaccine        = (bool) args[0];
    const double vac_efficacy = args[1];
    const double vac_coverage = args[2]; 
    const int vac_mech        = (int) args[3];
    const int vac_day         = (int) args[4];
    //const double vac_soc_dist = args[5]; // social distancing in effect post vaccination campaign

    Community* community = build_community(par);

    if (vaccine) {
        //double target_coverage  = coverage;
        double catchup_coverage = vac_coverage;
        const int target = 0;
        const int catchup_to = NUM_AGE_CLASSES - 1;

        for (int catchup_age = target; catchup_age <= catchup_to; catchup_age++) {
            //const int vacDate = julian_to_sim_day(par, JULIAN_TALLY_DATE + 1, RESTART_BURNIN);
            par->catchupVaccinationEvents.emplace_back(catchup_age, vac_day, catchup_coverage);
        }

//        par->vaccineTargetAge = target;
//        par->vaccineTargetCoverage = target_coverage;
//        par->vaccineTargetStartDate = julian_to_sim_day(par, JULIAN_TALLY_DATE + 1, RESTART_BURNIN);
//        par->seroTestFalsePos = seroTestFalsePos;
//        par->seroTestFalseNeg = seroTestFalseNeg;
    }

    if (vac_mech == 0) {
        par->VES                   = vac_efficacy;
        par->VES_NAIVE             = vac_efficacy;
        //par->VEH                    = 0.803; // fraction of hospitalized cases prevented by vaccine
        //par->linearlyWaningVaccine   = true;
        //par->vaccineImmunityDuration = 2*365;
        par->vaccineLeaky           = false;
        par->numVaccineDoses         = 1;
        //par->vaccineDoseInterval     = 182;
        //par->vaccineBoosting         = false;
        //par->vaccineSeroConstraint   = VACCINATE_ALL_SERO_STATUSES;
        //par->whoDiseaseOutcome       = INC_NUM_INFECTIONS;
    } else {
        cerr << "Unsupported vaccine mechanism: " << vac_mech << endl;
        exit(1);
    }

    seed_epidemic(par, community);
    simulate_epidemic(par, community, process_id);
    time (&end);
    double dif = difftime (end,start);

//    const int desired_intervention_output = FORECAST_DURATION - 1;
    vector<double> metrics = tally_counts(par, community, 0);

//    const vector<size_t> infections       = community->getNumNewlyInfected();
//    const vector<size_t> reported_cases   = community->getNumDetectedCasesReport();
//    const vector<size_t> hospitalizations = community->getNumHospPrev();
//    const vector<size_t> icu              = community->getNumIcuInc();
//    const vector<size_t> deaths           = community->getNumDetectedDeaths();
//    const vector<float> lockdown          = community->getTimedIntervention(NONESSENTIAL_BUSINESS_CLOSURE);

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

void usage() {
    cerr << "\n\tUsage: ./abc_sql abc_config_sql.json --process\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate -n <number of simulations per database write>\n\n";

}


int main(int argc, char* argv[]) {
    if (not (argc == 3 or argc == 5 or argc == 6) ) {
        usage();
        exit(100);
    }

    bool process_db = false;
    bool simulate_db = false;
    int buffer_size = -1;

    for (int i=2; i < argc;  i++ ) {
        if ( strcmp(argv[i], "--process") == 0  ) {
            process_db = true;
        } else if ( strcmp(argv[i], "--simulate") == 0  ) {
            simulate_db = true;
            buffer_size = buffer_size == -1 ? 1 : buffer_size;
        } else if ( strcmp(argv[i], "-n" ) == 0 ) {
            buffer_size = atoi(argv[++i]);
        } else {
            usage();
            exit(101);
        }
    }

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(string(argv[1]));
    if (process_db) {
        gsl_rng_set(RNG, time(NULL) * getpid()); // seed the rng using sys time and the process id
        abc->process_database(RNG);
    }

    if (simulate_db) {
        time(&GLOBAL_START_TIME);
        abc->set_simulator(simulator);
        abc->simulate_next_particles(buffer_size);
    }

    return 0;
}
