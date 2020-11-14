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
const string SIM_POP = "escambia"; // recent good T value == 0.042
//const string SIM_POP = "dade"; // recent good T value == 0.03
//const string SIM_POP = "florida";
const string HOME_DIR(std::getenv("HOME"));
const string pop_dir = HOME_DIR + "/work/covid-abm/pop/" + SIM_POP;
const string output_dir("/ufrc/longini/tjhladish/");
//const string imm_dir(output_dir + "");

const int RESTART_BURNIN          = 0;
const int FORECAST_DURATION       = 240;
const bool RUN_FORECAST           = true;
const int TOTAL_DURATION          = RUN_FORECAST ? RESTART_BURNIN + FORECAST_DURATION : RESTART_BURNIN;
const size_t JULIAN_TALLY_DATE    = 146; // intervention julian date - 1
const size_t JULIAN_START_YEAR    = 2020;
const double DEATH_UNDERREPORTING = 11807.0/20100.0; // FL Mar15-Sep5, https://www.nytimes.com/interactive/2020/05/05/us/coronavirus-death-toll-us.html

int to_sim_day(int julian_start, string date) { return Date::to_julian_day(date) - julian_start; }

//Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed) {
Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const string /*process_id*/) {
    Parameters* par = new Parameters();
    par->define_defaults();
    //const MmodsScenario ms = (MmodsScenario) args[0];
    const MmodsScenario ms = NUM_OF_MMODS_SCENARIOS; // if NUM_OF_MMODS_SCENARIOS, just run a normal simulation
    assert(ms <= NUM_OF_MMODS_SCENARIOS);
    if (ms < NUM_OF_MMODS_SCENARIOS) { par->numSurveilledPeople = 1e5; } // default is INT_MAX
    par->serial = serial;

    const float T = args[0]; //0.2;

    par->household_transmissibility   = T;
    par->social_transmissibility      = T*2.0; // assumes complete graphs between households
    par->workplace_transmissibility   = T;//10.0;
    par->school_transmissibility      = T;//5.0;
    par->hospital_transmissibility    = T;//10.0;
    par->nursinghome_transmissibility = T;
    par->randomseed              = rng_seed;
    par->dailyOutput             = false; // turn on for daily prevalence figure, probably uncomment filter in simulator.h for daily output to get only rel. days
    par->periodicOutput          = false;
    par->periodicOutputInterval  = 20;
    par->weeklyOutput            = false;
    par->monthlyOutput           = true;
    par->yearlyOutput            = true;
    par->abcVerbose              = false; // needs to be false to get WHO daily output
    par->julianYear              = JULIAN_START_YEAR;
    par->startDayOfYear          = Date::to_julian_day("2020-02-10");;
    par->runLength               = TOTAL_DURATION;
    par->annualIntroductionsCoef = 1;

    par->traceContacts = true; // needed for Rt calculation
    par->mmodsScenario = ms;      // MMODS_CLOSED, MMODS_2WEEKS, MMODS_1PERCENT, MMODS_OPEN

    {
        // Create model for how outcome-dependent detection probabilities change over time
        //par->reportedFraction = {0.0, 0.01, 0.5, 0.8, 1.0};      // fraction of asymptomatic, mild, severe, critical, and deaths reported
        //par->reportedFraction = {0.0, 0.2, 0.75, 0.75, 0.75};      // fraction of asymptomatic, mild, severe, critical, and deaths reported
        //par->probFirstDetection = {0.0, 0.12, 0.55, 0.1, 0.01};      // probability of being detected while {asymp, mild, severe, crit, dead} if not detected previously

        // probability of being detected while {asymp, mild, severe, crit, dead} if not detected previously
        const vector<double> initial_vals = {0.0, 0.1, 0.7, 0.1, 0.01};
        const vector<double> final_vals   = {0.05, 0.5, 0.3, 0.1, 0.01};
        const int isd = to_sim_day(par->startDayOfYear, "2020-06-01");
        const vector<int> inflection_sim_day = {isd, isd, isd, isd, isd};
        const vector<double> slopes = {0.1, 0.1, 0.1, 0.1, 0.1}; // sign is determined based on initial/final values

        par->createDetectionModel(initial_vals, final_vals, inflection_sim_day, slopes);

        vector<double> reported_frac_init  = par->toReportedFraction(initial_vals);
        vector<double> reported_frac_final = par->toReportedFraction(final_vals);
        cerr_vector(reported_frac_init); cerr << endl;
        cerr_vector(reported_frac_final); cerr << endl;

        //cerr << "Detection probability by day\n";
        //for (size_t d = 0; d < par->probFirstDetection.size(); ++d) {
        //    cerr << d; for (auto v: par->probFirstDetection[d]) cerr << " " << v; cerr << endl;
        //}
    }

    // These are only initial values for time-structured interventions.  They can be changed dynamically.
//    par->timedInterventions[SCHOOL_CLOSURE].resize(par->runLength, 0.0);
//    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(par->runLength, 0.0);
//    par->timedInterventions[SOCIAL_DISTANCING].resize(par->runLength, 0.0);

    par->timedInterventions[SCHOOL_CLOSURE].clear();
    par->timedInterventions[SCHOOL_CLOSURE].resize(to_sim_day(par->startDayOfYear, "2020-03-15"), 0.0);
    const size_t aug31 = to_sim_day(par->startDayOfYear, "2020-08-31");
    const size_t school_closed_duration = aug31 < par->runLength ? aug31 : par->runLength;
    par->timedInterventions[SCHOOL_CLOSURE].resize(school_closed_duration, 1.0);
    //par->timedInterventions[SCHOOL_CLOSURE].resize(par->runLength, 0.5); // 50% reopening on Aug 31
    par->timedInterventions[SCHOOL_CLOSURE].resize(par->runLength, 0.0); // 50% reopening on Aug 31

    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].clear();
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(to_sim_day(par->startDayOfYear, "2020-04-03"), 0.0);
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(to_sim_day(par->startDayOfYear, "2020-05-04"), 1.0);
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(par->runLength, 0.0);

    //const float mobility_logit_shift   = -2.5;
    //const float mobility_logit_stretch = 3.0;
    const float mobility_logit_shift   = 0.0;//-2.5;
    const float mobility_logit_stretch = -1.0;//3.0;
    // using createSocialDistancingModel() defines timedInterventions[SOCIAL_DISTANCING] values
    // NB: startDayOfYear, runLength, and julianYear must be defined in par before a social distancing model can be meaningfully created!
    // TODO - make that dependency something the user doesn't have to know or think about
    //par->createSocialDistancingModel(pop_dir + "/safegraph_mobility_index.csv", mobility_logit_shift, mobility_logit_stretch);
    // plot(as.Date(d$date), shiftstretch(1-d$smooth_ma7, shift=-2.7, stretch=2), ylim=c(0,1), type='l')
    //par->createSocialDistancingModel(pop_dir + "/../florida/mobility_comp-florida.csv", 1, mobility_logit_shift, mobility_logit_stretch);
    //par->createSocialDistancingModel(pop_dir + "/sgmi_escambia_comp.csv", mobility_logit_shift, mobility_logit_stretch);
    //par->createSocialDistancingModel(pop_dir + "/sgmi_dade_comp.csv", mobility_logit_shift, mobility_logit_stretch);
    //par->timedInterventions[SOCIAL_DISTANCING].resize(to_sim_day(par->startDayOfYear, "2020-05-15"));
    //par->timedInterventions[SOCIAL_DISTANCING].resize(par->runLength, 0.3);

    const size_t sd_metric_col_idx = (size_t) args[1];
    par->createSocialDistancingModel(pop_dir + "/../florida/gleam_metrics.csv", sd_metric_col_idx, mobility_logit_shift, mobility_logit_stretch);

    //par->defaultReportingLag = 14;
    par->createReportingLagModel(pop_dir + "/case_report_delay.csv");
    par->symptomToTestLag = 2;
    par->deathReportingLag = 9;

    const double max_icu_mortality_reduction = 0.5; // primarily due to use of dexamethasone
    const size_t icu_mortality_inflection_sim_day = to_sim_day(par->startDayOfYear, "2020-06-15");
    const double icu_mortality_reduction_slope = 0.5; // 0.5 -> change takes ~2 weeks; 0.1 -> ~2 months
    par->createIcuMortalityReductionModel(max_icu_mortality_reduction, icu_mortality_inflection_sim_day, icu_mortality_reduction_slope);
    par->icuMortalityFraction = 0.2;                        // to be fit; fraction of all deaths that occur in ICUs;
                                                            // used for interpreting empirical mortality data, *not within simulation*
    par->pathogenicityReduction = 0.6;                      // to be fit; fraction of infections missed in pathogenicity studies
                                                            // used for interpreting input data, *not within simulation*
//    par->susceptibilityCorrection = 1.0;
    par->define_susceptibility_and_pathogenicity();

    par->daysImmune = 730;
    par->VES = 0.0;

    //par->hospitalizedFraction = 0.25; // fraction of cases assumed to be hospitalized
    par->probInitialExposure = {3.0e-04};
    par->probDailyExposure   = {0.0e-05};

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


vector<double> tally_counts(const Parameters* par, Community* community, int pre_intervention_output) {
    const int discard_days = julian_to_sim_day(par, JULIAN_TALLY_DATE, RESTART_BURNIN-pre_intervention_output);

    //vector< vector<int> > severe      = community->getNumSevereCases();
    vector<size_t> symptomatic = community->getNumNewlySymptomatic();

    //vector< vector<int> > infected    = community->getNumNewlyInfected();
    const int num_years = FORECAST_DURATION + pre_intervention_output - 1; // typically 55
    vector<size_t> s_tally(num_years+1, 0);

    vector<double> metrics(num_years, 0.0);
    for (size_t t=discard_days; t<par->runLength; t++) {
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

vector<double> calc_Rt_moving_average(vector<pair<size_t, double>> Rt_pairs, size_t window) {
    assert(window % 2 == 1);
    vector<double> Rt_ma(Rt_pairs.size(), 0.0);

    size_t halfwindow = (window - 1)/2;
    for (size_t i = 0; i < Rt_ma.size(); ++i) {
        size_t window_ct = 0;
        double infections = 0;
        for (size_t wi = halfwindow >= i ? 0 : i - halfwindow;
             wi < (i + halfwindow + 1 >= Rt_pairs.size() ? Rt_pairs.size() : i + halfwindow + 1);
             ++wi) {
            const size_t ct = Rt_pairs[wi].first;
            const double Rt = Rt_pairs[wi].second;
            window_ct += ct;
            infections += ct*Rt;
        }
        Rt_ma[i] = infections / window_ct;
//        cerr << "day, inf, ct, Rt_ma, Rt_count, Rt_val: " << i << " " << infections << " " << window_ct << " " << Rt_ma[i] << " " << Rt_pairs[i].first << " " << Rt_pairs[i].second << endl;
    }
    return Rt_ma;
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

    Community* community = build_community(par);

    seed_epidemic(par, community);
    vector<string> plot_log_buffer = simulate_epidemic(par, community, process_id);
    vector<pair<size_t, double>> Rt = community->getMeanNumSecondaryInfections();
    vector<double> Rt_ma = calc_Rt_moving_average(Rt, 7);

    assert(Rt.size()+1 == plot_log_buffer.size()); // there's a header line
    for (size_t i = 1; i < plot_log_buffer.size(); ++i) {
        //plot_log_buffer[i] = plot_log_buffer[i] + "," + to_string(Rt[i-1].second);
        plot_log_buffer[i] = plot_log_buffer[i] + "," + to_string(Rt_ma[i-1]);
    }
    bool overwrite = true;
    write_daily_buffer(plot_log_buffer, process_id, "plot_log.csv", overwrite);
    int retval = system("Rscript simvis.R");
    if (retval == -1) { cerr << "System call to `Rscript simvis.R` failed\n"; }

    time (&end);
    double dif = difftime (end,start);

    map<string, vector<vector<size_t>>> ages_by_outcome;
    ages_by_outcome["cases"].resize(2*par->runLength/7);
    ages_by_outcome["hosp"].resize(2*par->runLength/7);
    ages_by_outcome["deaths"].resize(2*par->runLength/7);
    // mean and median age of (actual) cases, hospitalizations, and deaths over time
    vector<double> infections(NUM_AGE_CLASSES, 0.0);
    vector<double> cases(NUM_AGE_CLASSES, 0.0);
    vector<double> severe(NUM_AGE_CLASSES, 0.0);
    vector<double> critical(NUM_AGE_CLASSES, 0.0);
    vector<double> deaths(NUM_AGE_CLASSES, 0.0);

    for (Person* p: community->getPeople()) {
        if (p->getNumNaturalInfections() == 0) continue;
        const size_t age = p->getAge();
        infections[age]++;
        Infection* inf = p->getInfection();
//        if (inf->symptomatic()) cerr << "case week: " << inf->getSymptomTime()/7 << endl;
//        if (inf->hospital())    cerr << "hosp week: " << inf->getHospitalizedTime()/7 << endl;
//        if (inf->fatal())       cerr << "dead week: " << inf->getDeathTime()/7 << endl;
        if (inf->symptomatic()) { ages_by_outcome["cases"][inf->getSymptomTime()/7].push_back(age); cases[age]++; }
        if (inf->hospital())    { ages_by_outcome["hosp"][inf->getHospitalizedTime()/7].push_back(age); }
        if (inf->severe())      { severe[age]++; }
        if (inf->critical())    { critical[age]++; }
        if (inf->fatal())       { ages_by_outcome["deaths"][inf->getDeathTime()/7].push_back(age); deaths[age]++; }
    }
/*    cerr << "week case_mean case_median hosp_mean hosp_median death_mean death_median\n";
    for (size_t week = Date::julianWeek(par->startDayOfYear); week < par->runLength/7; ++week) {
        cerr << week << " "
             << mean(ages_by_outcome["cases"][week]) << " "
             << median(ages_by_outcome["cases"][week]) << " "
             << mean(ages_by_outcome["hosp"][week]) << " "
             << median(ages_by_outcome["hosp"][week]) << " "
             << mean(ages_by_outcome["deaths"][week]) << " "
             << median(ages_by_outcome["deaths"][week]) << endl;
    }*/

//    cerr << "age death critical severe symptomatic infection\n";
//    for (size_t age = 0; age < infections.size(); ++age) {
//            cerr << age << " " << deaths[age] << " " << critical[age] << " " << severe[age]  << " " << cases[age] << " " << infections[age]  << endl;
//    }

//    const int pre_intervention_output = 5; // years
//    const int desired_intervention_output = FORECAST_DURATION - 1;
    vector<double> metrics(0);// = tally_counts(par, community, pre_intervention_output);

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
    cerr << "\n\tUsage: ./abc_sql <MMOD_scenario> <rng_seed> <transmissibility>" << endl;
/*    cerr << "\n\tUsage: ./abc_sql abc_config_sql.json --process\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate -n <number of simulations per database write>\n\n";
*/
}


int main(int argc, char* argv[]) {
    if (argc != 4) { usage(); exit(-1); }
//    assert(argc == 4);
    vector<double> sim_args = {atof(argv[1]), atof(argv[2])};
    size_t seed = atol(argv[3]);
    size_t serial = 0;

    simulator(sim_args, seed, serial);

    return 0;
}
