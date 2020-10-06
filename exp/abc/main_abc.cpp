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
const string SIM_POP = "escambia"; // cycle through subpops when fitting?
//const string SIM_POP = "florida";
const string HOME_DIR(std::getenv("HOME"));
const string pop_dir = HOME_DIR + "/work/covid-abm/pop/" + SIM_POP;
const string output_dir("/ufrc/longini/tjhladish/");
//const string imm_dir(output_dir + "");

//const int RESTART_BURNIN       = 0;
//const int FORECAST_DURATION    = 364; // update to reflect amount of data we have to fit to
//const bool RUN_FORECAST        = true;
//const int TOTAL_DURATION       = RUN_FORECAST ? RESTART_BURNIN + FORECAST_DURATION : RESTART_BURNIN;
//const size_t JULIAN_TALLY_DATE = 146; // intervention julian date - 1
const size_t JULIAN_START_YEAR = 2020;
const size_t NUM_FITTED_WEEKS = 28;
const size_t FIRST_FITTED_JULIAN_DAY = Date::to_julian_day("2020-03-02");
const double DEATH_UNDERREPORTING = 11807.0/20100.0; // FL Mar15-Sep5, https://www.nytimes.com/interactive/2020/05/05/us/coronavirus-death-toll-us.html

//Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed) {
Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const string /*process_id*/) {
    Parameters* par = new Parameters();
    par->define_defaults();
    par->serial = serial;

    // ABC pars
    // 0: transmissibility
    // 1: start_date // julian day
    // 2: external_exposures
    // 3: first_detection_mild
    // 4: first_detection_severe
    // 5: icu_prob_given_death
    // 6: mobility_logit_shift
    // 7: mobility_logit_stretch

    const float T = args[0]; //0.25; // fit

    par->household_transmissibility   = T;
    par->social_transmissibility      = T/10.0; // assumes complete graphs between households
    par->workplace_transmissibility   = T;
    par->school_transmissibility      = T;
    par->hospital_transmissibility    = T/10.0;
    par->nursinghome_transmissibility = 2*T;
    //hospitalizedFraction = {0.0, 0.15, 0.9};
    //par->reportedFraction = {0.0, 0.01, 0.5, 0.8, 1.0};      // fraction of asymptomatic, mild, severe, critical, and deaths reported
    //par->reportedFraction = {0.0, 0.2, 0.75, 0.75, 0.75};      // fraction of asymptomatic, mild, severe, critical, and deaths reported
    //par->probFirstDetection = {0.0, 0.15, 0.6, 0.1, 0.01};      // probability of being detected while {asymp, mild, severe, crit, dead} if not detected previously
    par->probFirstDetection = {0.0, args[3], args[4], 0.1, 0.01};      // probability of being detected while {asymp, mild, severe, crit, dead} if not detected previously
    par->randomseed              = rng_seed;
    par->dailyOutput             = false; // turn on for daily prevalence figure, probably uncomment filter in simulator.h for daily output to get only rel. days
    par->periodicOutput          = false;
    par->periodicOutputInterval  = 20;
    par->weeklyOutput            = false;
    par->monthlyOutput           = true;
    par->yearlyOutput            = true;
    par->abcVerbose              = false; // needs to be false to get WHO daily output
    par->julianYear              = JULIAN_START_YEAR;
    par->startDayOfYear          = (int) args[1]; //45;
    par->runLength               = FIRST_FITTED_JULIAN_DAY - par->startDayOfYear + 7*NUM_FITTED_WEEKS + 1;
    par->annualIntroductionsCoef = 1;

    par->traceContacts = true;

    // These are only initial values for time-structured interventions.  They can be changed dynamically.
    par->timedInterventions[SCHOOL_CLOSURE].clear();
    par->timedInterventions[SCHOOL_CLOSURE].resize(30, 0.0);
    const size_t aug31 = 243 - par->startDayOfYear;
    const size_t school_closed_duration = aug31 < par->runLength ? aug31 : par->runLength;
    par->timedInterventions[SCHOOL_CLOSURE].resize(school_closed_duration, 1.0);
    par->timedInterventions[SCHOOL_CLOSURE].resize(par->runLength, 0.5); // 50% reopening on Aug 31

    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].clear();
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(45, 0.0);
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(75, 1.0);
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(par->runLength, 0.0);

    const float mobility_logit_shift   = args[6]; //-2.7;
    const float mobility_logit_stretch = args[7]; //2;
    // using createSocialDistancingModel() defines timedInterventions[SOCIAL_DISTANCING] values
    // NB: startDayOfYear, runLength, and julianYear must be defined in par before a social distancing model can be meaningfully created!
    // TODO - make that dependency something the user doesn't have to know or think about
    //par->createSocialDistancingModel(pop_dir + "/safegraph_mobility_index.csv", mobility_logit_shift, mobility_logit_stretch);
    // plot(as.Date(d$date), shiftstretch(1-d$smooth_ma7, shift=-2.7, stretch=2), ylim=c(0,1), type='l')
    par->createSocialDistancingModel(pop_dir + "/../florida/sgmi_fl_comp.csv", mobility_logit_shift, mobility_logit_stretch);
    //par->createSocialDistancingModel(pop_dir + "/sgmi_escambia_comp.csv", mobility_logit_shift, mobility_logit_stretch);
    //par->createSocialDistancingModel(pop_dir + "/sgmi_dade_comp.csv", mobility_logit_shift, mobility_logit_stretch);

    //par->defaultReportingLag = 14;
    par->createReportingLagModel(pop_dir + "/case_report_delay.csv");
    par->symptomToTestLag = 2;
    par->deathReportingLag = 7;

    const double max_icu_mortality_reduction = args[8]; // primarily due to use of dexamethasone
    //const size_t icu_mortality_inflection_sim_day = Date::to_julian_day("2020-06-25") - par->startDayOfYear;
    const size_t icu_mortality_inflection_sim_day = Date::to_julian_day("2020-07-01") - par->startDayOfYear;
    const double icu_mortality_reduction_slope = 0.1; // 0.5 -> change takes ~2 weeks; 0.1 -> ~2 months
    par->createIcuMortalityReductionModel(max_icu_mortality_reduction, icu_mortality_inflection_sim_day, icu_mortality_reduction_slope);
    par->icuMortalityFraction = args[5]; //0.5;             // to be fit; fraction of all deaths that occur in ICUs;
                                                            // used for interpreting empirical mortality data, *not within simulation*
    par->define_susceptibility_and_pathogenicity();

    par->daysImmune = 730;
    par->VES = 0.0;

    //par->hospitalizedFraction = 0.25; // fraction of cases assumed to be hospitalized
    par->probInitialExposure = {1.0e-04}; // for fit, assume exposures of ~1 per 50k people?
    par->probDailyExposure   = {args[2]}; //{1.0e-05};

    par->populationFilename       = pop_dir    + "/population-"         + SIM_POP + ".txt";
    par->comorbidityFilename      = pop_dir    + "/comorbidity-"        + SIM_POP + ".txt";
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


vector<double> tally_counts(const Parameters* par, Community* community) {
    // time series data currently start on March 2, 2020
    const size_t discard_days = FIRST_FITTED_JULIAN_DAY - par->startDayOfYear;
    //vector< vector<int> > severe      = community->getNumSevereCases();
    //vector<size_t> symptomatic = community->getNumNewlySymptomatic();
    const vector<size_t> all_reported_cases = community->getNumDetectedCasesReport();
    const vector<size_t> reported_deaths = community->getNumDetectedDeaths();


    //vector<size_t> infected    = community->getNumNewlyInfected();
    assert(par->runLength >= discard_days + NUM_FITTED_WEEKS*7);

    vector<double> metrics(2*NUM_FITTED_WEEKS, 0.0);
    for (size_t t=discard_days; t<par->runLength; t++) {
        const size_t w = (t-discard_days)/7;
        metrics[w] += all_reported_cases[t];
        metrics[NUM_FITTED_WEEKS + w] += reported_deaths[t]*DEATH_UNDERREPORTING;
    }

    // rescale metrics to be per 10k population
    for (auto& val: metrics) { val *= 1e4/community->getNumPeople(); }

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

    // ABC pars
    // 0: transmissibility
    // 1: start_date
    // 2: external_exposures
    // 3: first_detection_mild
    // 4: first_detection_severe
    // 5: icu_prob_given_death
    // 6: mobility_logit_shift
    // 7: mobility_logit_stretch

    Parameters* par = define_simulator_parameters(args, rng_seed, serial, process_id);
    vector<double> reported_frac = par->toReportedFraction(par->probFirstDetection);
    cerr_vector(reported_frac); cerr << endl;
    Community* community = build_community(par);

    seed_epidemic(par, community);
    vector<string> plot_log_buffer = simulate_epidemic(par, community, process_id);
    //vector<pair<size_t, double>> Rt = community->getMeanNumSecondaryInfections();
/*
    assert(Rt.size()+1 == plot_log_buffer.size()); // there's a header line
    for (size_t i = 1; i < plot_log_buffer.size(); ++i) {
        plot_log_buffer[i] = plot_log_buffer[i] + "," + to_string(Rt[i-1].second);
    }
    bool overwrite = true;
    write_daily_buffer(plot_log_buffer, process_id, "plot_log.csv", overwrite);
    int retval = system("Rscript simvis.R");
    if (retval == -1) { cerr << "System call to `Rscript simvis.R` failed\n"; }
*/
    time (&end);
    double dif = difftime (end,start);

//    const int desired_intervention_output = FORECAST_DURATION - 1;
    vector<double> metrics = tally_counts(par, community);

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
