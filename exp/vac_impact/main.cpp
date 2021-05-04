#include <unistd.h>
#include "AbcSmc.h"
#include "simulator.h"
#include <cstdlib>
#include "CCRC32.h"
#include "Utility.h"
#include <unordered_set>
#include <math.h>

using namespace std;

using covid::util::to_string;
using covid::util::mean;
using covid::util::stdev;
using covid::util::max_element;

time_t GLOBAL_START_TIME;

string calculate_process_id(vector<double> &args, string &argstring);
//const string SIM_POP = "escambia"; // recent good T value == 0.0253
//const string SIM_POP = "dade"; // recent good T value == 0.03
//const string SIM_POP = "florida";
//const string SIM_POP = "pseudo-1000k"; // recent good T value == 0.042
const string SIM_POP = "pseudo-300K"; // recent good T value == 0.042
const string HOME_DIR(std::getenv("HOME"));
const string pop_dir = HOME_DIR + "/work/covid-abm/pop/" + SIM_POP;
const string output_dir("/ufrc/longini/tjhladish/");
//const string imm_dir(output_dir + "");

const int RESTART_BURNIN          = 0;
const int FORECAST_DURATION       = 574;
//const int FORECAST_DURATION       = 330;
const int OVERRUN                 = 14; // to get accurate Rt estimates near the end of the forecast duration
const bool RUN_FORECAST           = true;
const int TOTAL_DURATION          = RUN_FORECAST ? RESTART_BURNIN + FORECAST_DURATION + OVERRUN : RESTART_BURNIN;
//const size_t JULIAN_TALLY_DATE    = 146; // intervention julian date - 1
const size_t JULIAN_START_YEAR    = 2020;
//const double DEATH_UNDERREPORTING = 11807.0/20100.0; // FL Mar15-Sep5, https://www.nytimes.com/interactive/2020/05/05/us/coronavirus-death-toll-us.html

//Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed) {
Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const string /*process_id*/) {
    Parameters* par = new Parameters();
    par->define_defaults();
    par->serial = serial;

    const float T = 0.055;//0.022; // 0.0215 for the fl pop, 0.022 for the pseudo 1000k pop

    par->household_transmissibility   = T;
    par->social_transmissibility      = T; // assumes complete graphs between households
    par->workplace_transmissibility   = T / 2.0;
    par->school_transmissibility      = T / 2.0;
    par->hospital_transmissibility    = T / 10.0;
    par->nursinghome_transmissibility = T * 3.0 / 2.0;
    par->randomseed              = rng_seed;
    par->dailyOutput             = false; // turn on for daily prevalence figure, probably uncomment filter in simulator.h for daily output to get only rel. days
    par->periodicOutput          = false;
    par->periodicOutputInterval  = 20;
    par->weeklyOutput            = false;
    par->monthlyOutput           = true;
    par->yearlyOutput            = true;
    par->abcVerbose              = false; // needs to be false to get WHO daily output
    par->startJulianYear         = JULIAN_START_YEAR;
    par->startDayOfYear          = Date::to_julian_day("2020-02-05");
    par->runLength               = TOTAL_DURATION;
    par->annualIntroductionsCoef = 1;

    const bool mutation          = (bool) args[5];

    vector<double> seasonality;
    for (size_t day = 0; day < 366; ++day) {
        // bi-modal curve corresponding to orlando, with seasonal forcing episilon = 0.25
        // need 366 days to accomodate leap years
        seasonality.push_back(1 - 0.25*sin((4*M_PI*(int) day)/366));
        //seasonality.push_back(1 - 0.15*2*(pow((1 + sin((4*M_PI*((int) day-60))/366))/2, 2)-0.5));
        //seasonality.push_back(1 - 0.25*2*(pow((1 + sin((4*M_PI*(day-60))/366))/2, 10.0)-0.5));
        //seasonality.push_back(1.0 - 0.25*((pow((1 + sin((4*M_PI*(day-60))/366))/2, 2) + (1 + sin((2*M_PI*(day+80))/366))/3.5)*(2/1.311193)-1));
        //seasonality.push_back(1.0);
    }
    cerr_vector(seasonality);
    par->seasonality = seasonality;

//    if (mutation) {
//        par->seasonality = vector<double>(Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2020-12-31"), 1.0);
//        for (int i = 0; i < 31; ++i ) { par->seasonality.push_back(1.0 + 0.5*i/31); }
//        par->seasonality.resize(par->runLength, 1.5);
//    } else {
//        par->seasonality.resize(par->runLength, 1.0);
//    }

    par->traceContacts = true; // needed for Rt calculation

    {
        // Create model for how outcome-dependent detection probabilities change over time
        //par->reportedFraction = {0.0, 0.01, 0.5, 0.8, 1.0};      // fraction of asymptomatic, mild, severe, critical, and deaths reported
        //par->reportedFraction = {0.0, 0.2, 0.75, 0.75, 0.75};      // fraction of asymptomatic, mild, severe, critical, and deaths reported
        //par->probFirstDetection = {0.0, 0.12, 0.55, 0.1, 0.01};      // probability of being detected while {asymp, mild, severe, crit, dead} if not detected previously

        const double RF_death_early = 0.8;//0.78; // overall probability of detecting death, at any point

        // probability of being detected while {asymp, mild, severe, crit, dead} if not detected previously
        vector<double> initial_vals = {0.01, 0.2, 0.7, 0.1};
        const double rho_death_ini  = 1.0 - (1.0 - RF_death_early)/((1.0 - initial_vals[0])*(1.0 - initial_vals[1])*(1.0 - initial_vals[2])*(1.0 - initial_vals[3]));
        initial_vals.push_back(rho_death_ini);

        vector<double> mid_vals   = {0.07, 0.6, 0.1, 0.1};
        const double rho_death_mid  = 1.0 - (1.0 - RF_death_early)/((1.0 - mid_vals[0])*(1.0 - mid_vals[1])*(1.0 - mid_vals[2])*(1.0 - mid_vals[3]));
        mid_vals.push_back(rho_death_mid);

        const double RF_death_late = 1.0;//0.78; // overall probability of detecting death, at any point
        vector<double> final_vals   = {0.3, 0.8, 0.0, 0.0};
        const double rho_death_fin  = 1.0 - (1.0 - RF_death_late)/((1.0 - final_vals[0])*(1.0 - final_vals[1])*(1.0 - final_vals[2])*(1.0 - final_vals[3]));
        final_vals.push_back(rho_death_fin);

        cerr << "init, mid, fin: " << rho_death_ini << " " << rho_death_mid << " " << rho_death_fin << endl;

        const int isd1 = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2020-06-01");
        const int isd2 = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2020-10-01");

        const vector<int> inflection1_sim_day = {isd1, isd1, isd1, isd1, isd1};
        const vector<int> inflection2_sim_day = {isd2, isd2, isd2, isd2, isd2};
        const vector<double> slopes = {0.1, 0.1, 0.1, 0.1, 0.1}; // sign is determined based on initial/final values

        par->createDetectionModel(initial_vals, mid_vals, final_vals, inflection1_sim_day, inflection2_sim_day, slopes, slopes);

        vector<double> reported_frac_init  = par->toReportedFraction(initial_vals);
        vector<double> reported_frac_mid   = par->toReportedFraction(mid_vals);
        vector<double> reported_frac_final = par->toReportedFraction(final_vals);
        cerr_vector(reported_frac_init); cerr << endl;
        cerr_vector(reported_frac_mid); cerr << endl;
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
    par->timedInterventions[SCHOOL_CLOSURE].resize(Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2020-03-15"), 0.0);
    const size_t aug31 = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2020-08-31");
    const size_t school_closed_duration = aug31 < par->runLength ? aug31 : par->runLength;
    par->timedInterventions[SCHOOL_CLOSURE].resize(school_closed_duration, 1.0);
    par->timedInterventions[SCHOOL_CLOSURE].resize(par->runLength, 0.5); // 50% reopening on Aug 31

    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].clear();
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2020-04-03"), 0.0);
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2020-05-04"), 1.0);
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(par->runLength, 0.0);

//    const float mobility_logit_shift   = -2.5;
//    const float mobility_logit_stretch = 2.0;
//    const float mobility_logit_shift   = 0.0;//-2.5;
//    const float mobility_logit_stretch = 1;//3.0;
    // using createSocialDistancingModel() defines timedInterventions[SOCIAL_DISTANCING] values
    // NB: startDayOfYear, runLength, and julianYear must be defined in par before a social distancing model can be meaningfully created!
    // TODO - make that dependency something the user doesn't have to know or think about
    //par->createSocialDistancingModel(pop_dir + "/safegraph_mobility_index.csv", mobility_logit_shift, mobility_logit_stretch);
    // plot(as.Date(d$date), shiftstretch(1-d$smooth_ma7, shift=-2.7, stretch=2), ylim=c(0,1), type='l')
    //par->createSocialDistancingModel(pop_dir + "/../florida/mobility_comp-florida.csv", 1, mobility_logit_shift, mobility_logit_stretch);
    //par->createSocialDistancingModel(pop_dir + "/sgmi_escambia_comp.csv", mobility_logit_shift, mobility_logit_stretch);
    //par->createSocialDistancingModel(pop_dir + "/sgmi_fl_comp.csv", 1, mobility_logit_shift, mobility_logit_stretch);
    //par->timedInterventions[SOCIAL_DISTANCING].resize(Date::to_sim_day(par->startDayOfYear, "2020-05-15"));
    //par->timedInterventions[SOCIAL_DISTANCING].resize(par->runLength, 0.3);

//    const size_t sd_metric_col_idx = 1;//(size_t) args[1];
//    par->createSocialDistancingModel(pop_dir + "/../gleam_metrics.csv", sd_metric_col_idx, mobility_logit_shift, mobility_logit_stretch);
//    for (auto &e: par->timedInterventions[SOCIAL_DISTANCING]) { e *= 0.8; }

    // Create social distancing model by linear interpolation ---
//    vector<TimeSeriesAnchorPoint> ap = { // tuning for whole FL model
//        {"2020-01-01", 0.0},
//        {"2020-03-10", 0.0},
//        {"2020-03-15", 0.1},
//        {"2020-04-01", 0.8},
//        {"2020-05-01", 0.7},
//        {"2020-06-01", 0.3},
//        {"2020-07-01", 0.2},
//        {"2020-08-01", 0.5},
//        {"2020-09-01", 0.5},
//        {"2020-11-15", 0.0}
//    };

    vector<TimeSeriesAnchorPoint> ap = { // tuning for whole FL model
        {"2020-01-01", 0.0},
//        {"2020-03-15", 0.1},
        {"2020-03-15", 0.0},
        {"2020-04-01", 0.6},
        {"2020-05-01", 0.6},
        {"2020-06-01", 0.2},
        {"2020-07-01", 0.3},
        {"2020-08-01", 0.3},
        {"2020-09-01", 0.3},
        {"2020-11-15", 0.2}
    };

    par->timedInterventions[SOCIAL_DISTANCING].clear();
  const string sim_start_date = Date::to_ymd(par->startJulianYear, par->startDayOfYear);

    par->timedInterventions[SOCIAL_DISTANCING] = Date::linInterpolateTimeSeries(ap, par->startJulianYear, par->startDayOfYear);
  //par->timedInterventions[SOCIAL_DISTANCING].insert(par->timedInterventions[SOCIAL_DISTANCING].begin(), v.begin(), v.end());

    const double last_value = par->timedInterventions[SOCIAL_DISTANCING].back();
    par->timedInterventions[SOCIAL_DISTANCING].resize(par->runLength, last_value);
    // -------------------------------------------------------------

    //par->defaultReportingLag = 14;
    par->createReportingLagModel(pop_dir + "/../case_report_delay.csv");
    par->symptomToTestLag = 2;
//    par->meanDeathReportingLag = 9;

    const double max_icu_mortality_reduction = 0.4;         // primarily due to use of dexamethasone
    const size_t icu_mortality_inflection_sim_day = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2020-06-15");
    const double icu_mortality_reduction_slope = 0.5;       // 0.5 -> change takes ~2 weeks; 0.1 -> ~2 months
    par->createIcuMortalityReductionModel(max_icu_mortality_reduction, icu_mortality_inflection_sim_day, icu_mortality_reduction_slope);
    par->icuMortalityFraction = 0.2;                        // to be fit; fraction of all deaths that occur in ICUs;
                                                            // used for interpreting empirical mortality data, *not within simulation*
    par->pathogenicityReduction = 0.65;                      // to be fit; fraction of infections missed in pathogenicity studies
                                                            // used for interpreting input data, *not within simulation*
//    par->susceptibilityCorrection = 1.0;
    par->define_susceptibility_and_pathogenicity();

    //par->daysImmune = 730; // changing this to be a function call
    par->VES = {0.0};

    //par->hospitalizedFraction = 0.25; // fraction of cases assumed to be hospitalized
    par->probInitialExposure = {5.0e-04};
    par->probDailyExposure   = {1.0e-05};

    par->populationFilename       = pop_dir    + "/population-"         + SIM_POP + ".txt";
    par->comorbidityFilename      = pop_dir    + "/comorbidity-"        + SIM_POP + ".txt";
    par->locationFilename         = pop_dir    + "/locations-"          + SIM_POP + ".txt";
    par->networkFilename          = pop_dir    + "/network-"            + SIM_POP + ".txt";
    par->publicActivityFilename   = pop_dir    + "/public-activity-"    + SIM_POP + ".txt";

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
    const size_t num_weeks = (par->runLength - discard_days - OVERRUN)/7;

    //vector<size_t> infected    = community->getNumNewlyInfected();
    //vector< vector<int> > severe      = community->getNumSevereCases();
    vector<size_t> symptomatic = community->getNumNewlySymptomatic();
    vector<size_t> dead        = community->getNumNewlyDead();

    // pair of num of primary infections starting this day, and mean num secondary infections they cause
    vector<pair<size_t, double>> R = community->getMeanNumSecondaryInfections();
    vector<size_t> Rt_incidence_tally(num_weeks, 0);

    vector<double> metrics(num_weeks*3, 0.0);
    for (size_t t = discard_days; t < par->runLength - OVERRUN; t++) {
        const size_t w = (t-discard_days)/7; // which reporting week are we in?
        metrics[w]                 += symptomatic[t];
        metrics[num_weeks + w]     += dead[t];
        metrics[2 * num_weeks + w] += R[t].first > 0 ? R[t].first*R[t].second : 0;
        Rt_incidence_tally[w]      += R[t].first;
    }

    for (size_t w = 0; w < num_weeks; ++w) {
        metrics[2 * num_weeks + w] /= Rt_incidence_tally[w] > 0 ? Rt_incidence_tally[w] : 1.0;
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
    cerr << "rng seed: " << rng_seed << endl;
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
// 0  "vaccine", "short_name" : "vac", "num_type"   : "INT", "par1"       : 0, "par2"       : 1, "step"       : 1}, 
// 1  "vaccine_efficacy_susceptibility", "short_name" : "VES", "num_type"   : "FLOAT", "par1"       : 0.6, "par2"       : 0.6, "step"       : 0.0}, 
// 2  "vaccine_efficacy_pathogenicity", "short_name" : "VEP", "num_type"   : "FLOAT", "par1"       : 0.95, "par2"       : 0.95, "step"       : 0.0}, 
// 3  "vaccination_rate", "short_name" : "vac_rate", "num_type"   : "FLOAT", "par1"       : 0.0012, "par2"       : 0.003, "step"       : 0.0018}, 
// 4  "realization", "num_type"   : "INT", "par1"       : 0, "par2"       : 99}, 
// 5  "mutation", "num_type"   : "INT", "par1"       : 0, "par2"       : 1, "step"       : 1}
    
    enum VacRate { SLOW, FAST, NUM_OF_VAC_RATES };

    Parameters* par = define_simulator_parameters(args, rng_seed, serial, process_id);
    const bool vaccine             = (bool) args[0];
    const double vac_efficacy_susc = args[1];
    const double vac_efficacy_path = args[2];
    const VacRate vac_rate         = (VacRate) args[3]; // 0 == slow, 1 == fast
    // const size_t realization    = (size_t) args[4];
    // const bool mutation         = (bool) args[5];

    Community* community = build_community(par);

    if (vaccine) {
        //double target_coverage  = coverage;
        double catchup_coverage = 0.5;
        const int target = 16;
        const int senior_threshold = 65;
        size_t senior_vac_campaign_duration = vac_rate == SLOW ? 86 : 35; // 0.12% per day vs 0.3% per day for [65, max_age]
        size_t general_vac_campaign_duration = vac_rate == SLOW ? 257 : 103; // same, for those [16, 64]
        const size_t senior_vac_sim_day  = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2021-01-11"); // jan 1 + 10 days for efficacy
        const size_t general_vac_sim_day = senior_vac_sim_day + senior_vac_campaign_duration;

        for (int catchup_age = senior_threshold; catchup_age <= NUM_AGE_CLASSES - 1; catchup_age++) {
        cerr << "scheduling vaccinations:" << catchup_age << " " << senior_vac_sim_day << endl;
            //const int vacDate = julian_to_sim_day(par, JULIAN_TALLY_DATE + 1, RESTART_BURNIN);
            par->catchupVaccinationEvents.emplace_back(senior_vac_sim_day, senior_vac_campaign_duration, catchup_age, catchup_coverage);
        }

        for (int catchup_age = target; catchup_age < senior_threshold; catchup_age++) {
        cerr << "scheduling vaccinations:" << catchup_age << " " << general_vac_sim_day << endl;
            //const int vacDate = julian_to_sim_day(par, JULIAN_TALLY_DATE + 1, RESTART_BURNIN);
            par->catchupVaccinationEvents.emplace_back(general_vac_sim_day, general_vac_campaign_duration, catchup_age, catchup_coverage);
        }

        par->VES                   = {vac_efficacy_susc};
        par->VES_NAIVE             = {vac_efficacy_susc};
        par->VEP                   = {vac_efficacy_path};
        par->vaccineLeaky          = false;
        par->numVaccineDoses       = 1;
    }

    seed_epidemic(par, community);
    vector<string> plot_log_buffer = simulate_epidemic(par, community, process_id);

// comment out this block if simvis.R is not needed
{
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
}

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
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate --serial <serial to run>\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate --posterior <index to run>\n\n";

}


int main(int argc, char* argv[]) {
    if (not (argc == 3 or argc == 5 or argc == 6) ) {
        usage();
        exit(100);
    }

    bool process_db = false;
    bool simulate_db = false;
    int buffer_size = -1;
    int requested_serial = -1;
    int requested_posterior_idx = -1;

    for (int i=2; i < argc;  i++ ) {
        if ( strcmp(argv[i], "--process") == 0  ) {
            process_db = true;
        } else if ( strcmp(argv[i], "--simulate") == 0  ) {
            simulate_db = true;
            buffer_size = buffer_size == -1 ? 1 : buffer_size;
        } else if ( strcmp(argv[i], "-n" ) == 0 ) {
            buffer_size = atoi(argv[++i]);
        } else if ( strcmp(argv[i], "--serial" ) == 0 ) {
            requested_serial = atoi(argv[++i]);
        } else if ( strcmp(argv[i], "--posterior" ) == 0 ) {
            requested_posterior_idx = atoi(argv[++i]);
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
        if (requested_serial > -1) {
            abc->simulate_particle_by_serial(requested_serial);
        } else if (requested_posterior_idx > -1) {
            abc->simulate_particle_by_posterior_idx(requested_posterior_idx);
        } else {
            abc->simulate_next_particles(buffer_size);
        }
    }

    return 0;
}
