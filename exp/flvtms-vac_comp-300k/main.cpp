#include <chrono>
#include <unistd.h>
#include "AbcSmc.h"
#include "simulator.h"
#include <cstdlib>
#include "CCRC32.h"
#include "Utility.h"
#include <math.h>
#include "../lib/exp_util.h"

#if __has_include("local.h")
#include "local.h"
#endif

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
//const string SIM_POP = "sim_pop-florida"; // recent good T value == 0.042
//const string HOME_DIR(std::getenv("HOME")); // commented out to test local path header

#ifdef LOCAL_HEADER
const std::string HOME_DIR = getHOME();
#else
const std::string HOME_DIR = std::getenv("HOME");
#endif

const string pop_dir = HOME_DIR + "/work/covid-abm/pop/" + SIM_POP;
const string output_dir("/ufrc/longini/tjhladish/");
//const string imm_dir(output_dir + "");
//const string vaccination_file = pop_dir + "/../fl_vac/fl_vac_v4.txt";

const int RESTART_BURNIN          = 0;
const int FORECAST_DURATION       = 747; // stop after omicron
//const int FORECAST_DURATION       = 468; // stop prior to delta
const int OVERRUN                 = 14; // to get accurate Rt estimates near the end of the forecast duration
const bool RUN_FORECAST           = true;
int TOTAL_DURATION          = RUN_FORECAST ? RESTART_BURNIN + FORECAST_DURATION + OVERRUN : RESTART_BURNIN;
//const size_t JULIAN_TALLY_DATE    = 146; // intervention julian date - 1
const size_t JULIAN_START_YEAR    = 2020;
//const double DEATH_UNDERREPORTING = 11807.0/20100.0; // FL Mar15-Sep5, https://www.nytimes.com/interactive/2020/05/05/us/coronavirus-death-toll-us.html
bool autotune                     = false;
const int FL_POP                  = 21538187;   // as of 2020 census

Parameters* define_simulator_parameters(vector<double> /*args*/, const unsigned long int rng_seed, const unsigned long int serial, const string /*process_id*/) {
    Parameters* par = new Parameters();
    par->define_defaults();
    par->serial = serial;

    const float T = 0.07; //0.022; // 0.0215 for the fl pop, 0.022 for the pseudo 1000k pop

    par->household_transmission_haz_mult   = T;
    par->social_transmission_haz_mult      = T / 2.0; // assumes complete graphs between households
    par->workplace_transmission_haz_mult   = T / 2.0;
    par->school_transmission_haz_mult      = T / 2.0;
    par->hospital_transmission_haz_mult    = T / 10.0;
    par->nursinghome_transmission_haz_mult = T * 3.0 / 2.0;
    par->randomseed              = rng_seed;
    par->dailyOutput             = false; // turn on for daily prevalence figure, probably uncomment filter in simulator.h for daily output to get only rel. days
    par->periodicOutput          = false;
    par->periodicOutputInterval  = 20;
    par->weeklyOutput            = false;
    par->monthlyOutput           = true;
    par->yearlyOutput            = true;
    par->abcVerbose              = false; // needs to be false to get WHO daily output
    par->startJulianYear         = JULIAN_START_YEAR;
    par->startDayOfYear          = Date::to_julian_day("2020-02-10");
    par->runLength               = TOTAL_DURATION;
    //par->annualIntroductionsCoef = 1;

    par->beginContactTracing           = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2021-05-01");
    par->contactTracingCoverage        = 0.7;
    par->contactTracingEV[HOME]        = 5.0;
    par->contactTracingEV[WORK]        = 3.0;
    par->contactTracingEV[SCHOOL]      = 3.0;
    par->contactTracingEV[HOSPITAL]    = 0.0;
    par->contactTracingEV[NURSINGHOME] = 5.0;
    par->contactTracingDepth           = 2;

    par->quarantineProbability  = {0.0, 0.0, 0.0};
    // par->quarantineProbability  = {0.9, 0.75, 0.5};
    par->quarantineDuration = 10;

    vector<double> seasonality;
    for (size_t day = 0; day < 366; ++day) {
        seasonality.push_back(1 - 0.15*sin((4*M_PI*((int) day-60))/366));
    }
    par->seasonality = seasonality;

    par->traceContacts = true; // needed for Rt calculation

    {
        // Create model for how outcome-dependent detection probabilities change over time
        const double RF_death = 1.0; // overall probability of detecting death, at any point

        // probability of being detected while {asymp, mild, severe, crit, dead} if not detected previously
        vector<vector<double>> first_detection_probs = { {0.0, 0.1, 0.7, 0.1}, // up to "2020-06-01"
                                                         {0.1, 0.5, 0.5, 0.1}, // up to "2020-10-01"
                                                         {0.3, 0.9, 0.5, 0.1}, // up to "2021-02-15"
                                                         {0.1, 0.7, 0.75, 0.1},// up to "2021-10-01"
                                                         {0.1, 0.6, 0.9, 0.1} };

        add_death_probabilities(first_detection_probs, RF_death);
        vector<string> inflection_dates = {"2020-06-01",
                                           "2020-10-01",
                                           "2021-02-15",
                                           "2021-10-01"};
        vector<vector<int>> inflection_matrix = create_sim_day_matrix(par, inflection_dates);
        // sign of slope is determined based on initial/final values
        vector<vector<double>> slope_matrix   = create_slope_matrix(inflection_dates.size(), 0.1); // assume same logistic slope for all reporting transitions

        par->createDetectionModel(first_detection_probs, inflection_matrix, slope_matrix);

        vector<vector<double>> RF_matrix = as_reported_fractions(par, first_detection_probs);
        cerr_matrix(RF_matrix);

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

    const size_t aug31_2020 = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2020-08-31");
    const size_t school_closed_duration = min(aug31_2020, par->runLength);
    par->timedInterventions[SCHOOL_CLOSURE].resize(school_closed_duration, 1.0);

    const size_t jun16_2021 = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2021-06-16");
    const size_t end_school_year_2020_2021 = min(jun16_2021, par->runLength);
    par->timedInterventions[SCHOOL_CLOSURE].resize(end_school_year_2020_2021, 0.5); // 50% reopening on Aug 31

    const size_t aug09_2021 = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2021-08-09");
    const size_t end_summer_2021 = min(aug09_2021, par->runLength);
    par->timedInterventions[SCHOOL_CLOSURE].resize(end_summer_2021, 1.0); //
    par->timedInterventions[SCHOOL_CLOSURE].resize(par->runLength, 0.2); // reopen after beinning of 2021-2022 school year

    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].clear();
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2020-04-03"), 0.0);
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2020-05-04"), 1.0);
    par->timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE].resize(par->runLength, 0.0);

/*    vector<TimeSeriesAnchorPoint> ap = { // tuning for whole FL model
        {"2020-01-01", 0.0},
        {"2020-03-10", 0.20},
        {"2020-03-15", 0.8},
        {"2020-04-01", 0.7},
        {"2020-05-01", 0.6},
        {"2020-06-01", 0.05},
        {"2020-07-01", 0.05},
        {"2020-08-01", 0.6},
        {"2020-09-01", 0.5},
        {"2020-10-01", 0.1},
        {"2020-11-01", 0.0},
        {"2020-12-01", 0.1},
        {"2021-01-01", 0.1},
        {"2021-02-01", 0.3},
        {"2021-03-01", 0.4},
        {"2021-04-01", 0.0},
        {"2021-05-01", 0.3},
        {"2021-06-01", 0.2},
        {"2021-07-01", 0.0},
        {"2021-08-01", 0.1},
        {"2021-09-01", 0.25},
        {"2021-10-01", 0.125},
        {"2021-11-01", 0.00}
    };

    par->timedInterventions[SOCIAL_DISTANCING].clear();
  const string sim_start_date = Date::to_ymd(par->startJulianYear, par->startDayOfYear);

    par->timedInterventions[SOCIAL_DISTANCING] = Date::linInterpolateTimeSeries(ap, par->startJulianYear, par->startDayOfYear);
  //par->timedInterventions[SOCIAL_DISTANCING].insert(par->timedInterventions[SOCIAL_DISTANCING].begin(), v.begin(), v.end());

    const double last_value = par->timedInterventions[SOCIAL_DISTANCING].back();
    par->timedInterventions[SOCIAL_DISTANCING].resize(par->runLength, last_value);
    */
    // -------------------------------------------------------------

    //par->defaultReportingLag = 14;
    par->createReportingLagModel(pop_dir + "/../case_report_delay.csv");
    par->symptomToTestLag = 2;
//    par->meanDeathReportingLag = 9;

    const double max_icu_mortality_reduction = 0.35;         // primarily due to use of dexamethasone
    const size_t icu_mortality_inflection_sim_day = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2020-06-15");
    const double icu_mortality_reduction_slope = 1.0;       // 0.5 -> change takes ~2 weeks; 0.1 -> ~2 months
    par->createIcuMortalityReductionModel(max_icu_mortality_reduction, icu_mortality_inflection_sim_day, icu_mortality_reduction_slope);
    par->icuMortalityFraction = 0.43;                       // 0.72*0.6; fraction of all deaths that occur in ICUs
                                                            // used for interpreting empirical mortality data, *not within simulation*
                                                            // 0.7229464 = fraction of covid FL deaths that were inpatient (1/2020 through 1/2022)
                                                            // --> https://www.cdc.gov/nchs/nvss/vsrr/covid_weekly/index.htm#PlaceDeath
                                                            // ~0.6 = fraction of inpatient deaths in US that were ventilated (proxy for ICU)
                                                            // --> https://www.cdc.gov/nchs/covid19/nhcs/hospital-mortality-by-week.htm

    par->pathogenicityReduction = 0.0;                      // to be fit; fraction of infections missed in pathogenicity studies
                                                            // used for interpreting input data, *not within simulation*
//    par->susceptibilityCorrection = 1.0;
    par->define_susceptibility_and_pathogenicity();

    //par->daysImmune = 730; // changing this to be a function call
    par->VES = {{WILDTYPE, {0.0}}};

    //par->hospitalizedFraction = 0.25; // fraction of cases assumed to be hospitalized
    par->probInitialExposure = {1.0e-04};
    //par->probDailyExposure   = vector(120, 2.0e-05);        // introductions are initially lower
    //par->probDailyExposure.resize(par->runLength, 2.0e-04); // and then pick up after ~ 4 months
    par->probDailyExposure   = {1.0e-04};        // introductions are initially lower

    par->populationFilename       = pop_dir    + "/population-"         + SIM_POP + ".txt";
    par->comorbidityFilename      = pop_dir    + "/comorbidity-"        + SIM_POP + ".txt";
    par->locationFilename         = pop_dir    + "/locations-"          + SIM_POP + ".txt";
    par->networkFilename          = pop_dir    + "/network-"            + SIM_POP + ".txt";
    par->publicActivityFilename   = pop_dir    + "/public-activity-"    + SIM_POP + ".txt";
    par->rCaseDeathFilename       = "./rcasedeath-florida.csv";
    par->vaccinationFilename      = "./state_based_counterfactual_doses.txt"; //"./counterfactual_doses_v2.txt";
    // par->doseFilename            = "./counterfactual_doses.txt"; //"./dose_data/FL_doses.txt"; //pop_dir    + "/../fl_vac/doses.txt";
    par->riskGroupsFilename       = "./300K_sample_pop_groups.txt";

    par->behavioral_autotuning = autotune;
    par->behavior_fitting_data_target = CASES;
    par->death_tuning_offset = 18; // 18 is median lag b/n infection and death; 8 is median lag b/n detection and death
    par->tuning_window = 14;
    par->num_preview_windows = 3;
    par->behaviorFilename = "autotuning_dataset_dump.csv";
    par->autotuningFilename = "autotuning_dataset_dump.csv";

    par->dump_simulation_data = false;

    return par;
}


void define_strain_parameters(Parameters* par) {
    par->IEP                   = 0.75;
    par->IEH                   = 0.5;

    par->VES                   = {{WILDTYPE, {0.40, 0.80}}, {ALPHA, {0.40, 0.80}}, {DELTA, {0.40, 0.80}}, {OMICRON, {0.40, 0.80}}}; // efficacy currently is being reduced in Person.cpp
    par->VEP                   = {{WILDTYPE, {0.67, 0.75}}, {ALPHA, {0.67, 0.75}}, {DELTA, {0.67, 0.75}}, {OMICRON, {0.67, 0.75}}};
    par->VEH                   = {{WILDTYPE, {0.9,  1.0}},  {ALPHA, {0.9,  1.0}},  {DELTA, {0.9,  0.93}}, {OMICRON, {0.48, 0.96}}};
    par->VEI                   = {{WILDTYPE, {0.4,  0.8}},  {ALPHA, {0.4,  0.8}},  {DELTA, {0.4,  0.8}},  {OMICRON, {0.2,  0.4}}};
    par->VEF                   = {{WILDTYPE, {0.0,  0.0}},  {ALPHA, {0.0,  0.0}},  {DELTA, {0.0,  0.0}},  {OMICRON, {0.0,  0.0}}}; // effect likely captured by VEH
    par->VES_NAIVE             = par->VES;

    par->strainPars[ALPHA].relInfectiousness   = 1.6;
    par->strainPars[ALPHA].relPathogenicity    = 1.1;
    par->strainPars[ALPHA].immuneEscapeProb    = 0.15;

    par->strainPars[DELTA].relInfectiousness   = par->strainPars[ALPHA].relInfectiousness * 1.6;
    par->strainPars[DELTA].relPathogenicity    = par->strainPars[ALPHA].relPathogenicity * 2.83;
    par->strainPars[DELTA].relSeverity         = 1.3; // relSeverity only applies if not vaccine protected; CABP - may be more like 1.3 based on mortality increase
    par->strainPars[DELTA].relIcuMortality     = 4.0; // TODO - this is due to icu crowding.  should be represented differently
    par->strainPars[DELTA].immuneEscapeProb    = 0.2;

    const double appxNonOmicronInfPd     = 9.0;
    const double appxOmicronInfPd        = 6.0;
    const double relInfectiousnessDenom  = (1.0 - pow(1.0 - par->household_transmission_haz_mult, appxOmicronInfPd/appxNonOmicronInfPd)) / par->household_transmission_haz_mult;

    par->strainPars[OMICRON].immuneEscapeProb  = 0.5;
    par->strainPars[OMICRON].relInfectiousness = par->strainPars[DELTA].relInfectiousness * 2.0 / relInfectiousnessDenom; // CABP: 2.148 would be justified
    par->strainPars[OMICRON].relPathogenicity  = par->strainPars[DELTA].relPathogenicity * 0.5;
    par->strainPars[OMICRON].relSeverity       = par->strainPars[DELTA].relSeverity * 0.75;
    par->strainPars[OMICRON].relIcuMortality   = 2.0;
    par->strainPars[OMICRON].symptomaticInfectiousPeriod = appxNonOmicronInfPd - 1;
    par->strainPars[OMICRON].relSymptomOnset = 0.5;     // roughly based on MMWR Early Release Vol. 70 12/28/2021

    cerr << "delta, omicron rel infectiousness: " << par->strainPars[DELTA].relInfectiousness << " " << par->strainPars[OMICRON].relInfectiousness << endl;

    par->crossProtectionMatrix = {{1, 0, 0, 0},    // WILDTYPE
                                  {0, 1, 0, 0},    // ALPHA
                                  {0, 0, 1, 0},    // DELTA
                                  {0, 0, 0, 1}};   // OMICRON
}


// REFACTOR parseVaccineFile()
void parseVaccineFile(const Parameters* par, Community* community, Vac_Campaign* vc, const size_t counterfactual_scenario, vector<bool> dose_pooling_flags,
                      bool adjust_std_to_bin_pop, bool adjust_urg_to_bin_pop) {
    // parses the JSON argument into a form to use in vax file parsing
    string counterfactual_reference_loc;
    const vector<string> loc_lookup = {"FL", "VT", "MS"};
    counterfactual_reference_loc = loc_lookup[counterfactual_scenario];

    // check that vaccinationFilename exists and can be opened
    ifstream iss(par->vaccinationFilename);
    string buffer;
    istringstream line;

    if (!iss) {
        cerr << "ERROR: vaccination file " << par->vaccinationFilename << " not found." << endl;
        exit(-1);
    }

    // variables for queue file processing (date ref_location bin_min bin_max dose n_doses_p10k)
    string date, ref_loc;
    int bin_min, bin_max, dose, is_urg;
    double doses_p10k;
    size_t first_vac_day = par->runLength;

    // save unique age bin boundaries encountered
    set<int> unique_bin_mins, unique_bin_maxs;

    // temmporary daily doses per 10k availability indexed by [day][dose][age bin]
    vector< vector< map<int, double> > > std_doses_in(par->runLength, vector< map<int, double> >(par->numVaccineDoses));
    vector< vector< map<int, double> > > urg_doses_in(par->runLength, vector< map<int, double> >(par->numVaccineDoses));

    // parse input file
    while (getline(iss, buffer)) {
        line.clear();
        line.str(buffer);

        if (line >> date >> ref_loc >> bin_min >> bin_max >> dose >> is_urg >> doses_p10k) {
            const size_t sim_day = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, date);
            if (sim_day < first_vac_day) { first_vac_day = sim_day; }

            // skip lines of data not pertaining to this counterfactual_reference_loc or are beyond runLength
            if (not ((ref_loc == counterfactual_reference_loc) and (sim_day < par->runLength))) { continue; }

            // will only insert values not encountered before
            unique_bin_mins.insert(bin_min);
            unique_bin_maxs.insert(bin_max);

            // save pop adjusted number of doses
            if (is_urg) {
                urg_doses_in[sim_day][dose - 1][bin_min] = doses_p10k;
            } else {
                std_doses_in[sim_day][dose - 1][bin_min] = doses_p10k;
            }
        }
    }
    iss.close();

    // TODO: do dose availability projections if runLength > last day of dose data file

    // set start of general campaign to first day of dose data
    vc->set_start_of_campaign(GENERAL_CAMPAIGN, first_vac_day);

    // call Vac_Campaign function to create age_bin_lookup data structure
    vc->generate_age_bins(community, unique_bin_mins, unique_bin_maxs);

    // for any day, dose, bin combinations with no data, fill with 0
    // otherwise, mulitply by proper age bin population
    map<int, int> bin_pops = vc->get_unique_age_bin_pops();
    vector< vector< map<int, int> > > std_doses_available(par->runLength, vector< map<int, int> >(par->numVaccineDoses));
    vector< vector< map<int, int> > > urg_doses_available(par->runLength, vector< map<int, int> >(par->numVaccineDoses));

    for (int day = 0; day < (int) par->runLength; ++day) {
        for (int bin : vc->get_unique_age_bins()) {
            for (int dose = 0; dose < par->numVaccineDoses; ++dose) {
                const double std_adj = adjust_std_to_bin_pop ? (double)bin_pops[bin]/1e4 : (double)community->getNumPeople()/1e4;
                const double urg_adj = adjust_urg_to_bin_pop ? (double)bin_pops[bin]/1e4 : (double)community->getNumPeople()/1e4;
                if (std_doses_in[day][dose].count(bin)) {
                    std_doses_available[day][dose][bin] = (int)(std_doses_in[day][dose][bin] * std_adj);
                } else {
                    std_doses_available[day][dose][bin] = 0;
                }

                if (urg_doses_in[day][dose].count(bin)) {
                    urg_doses_available[day][dose][bin] = (int)(urg_doses_in[day][dose][bin] * urg_adj);
                } else {
                    urg_doses_available[day][dose][bin] = 0;
                }
            }
        }
    }

    // save bin pop adjusted doses available to the vac_campaign
    vc->set_pool_urg_doses(dose_pooling_flags[0]);
    vc->set_pool_std_doses(dose_pooling_flags[1]);
    vc->set_pool_all_doses(dose_pooling_flags[2]);

    vc->init_doses_available(urg_doses_available, std_doses_available);
    // vc->set_doses_available(doses_available);
    // vc->set_urg_doses_available(urg_doses_available);
    // vc->init_orig_doses_available();
}

// REFACTOR generateVac_Campaign()
Vac_Campaign* generateVac_Campaign(const Parameters* par, Community* community, const size_t counterfactual_scenario, vector<bool> dose_pooling_flags,
                                   bool adjust_std_to_bin_pop, bool adjust_urg_to_bin_pop) {
    // create a new Vac_Campaign
    Vac_Campaign* vc = new Vac_Campaign(par);
    vc->set_rng(VAX_RNG);

    // parse input file to set daily doses available and generate unique, mutually exclusive age bins for the entire population
    cerr << "Reading vaccinations ... ";
    parseVaccineFile(par, community, vc, counterfactual_scenario, dose_pooling_flags, adjust_std_to_bin_pop, adjust_urg_to_bin_pop);
    cerr << "done." << endl;

    // initialiaze eligibility queue
    vc->init_eligibility_queue(community);

    // set newly initialized vac_campaign to community
    community->setVac_Campaign(vc);

    return vc;
}

int julian_to_sim_day (const Parameters* par, const size_t julian, const int intervention_year) {
    int startDate = intervention_year*365 + julian - par->startDayOfYear;
    if (julian < par->startDayOfYear) { // start intervention in following year
        startDate += 365;
    }
    return startDate;
}


vector<double> tally_counts(const Parameters* par, Community* community, int discard_days) {
    assert((int) par->runLength >= discard_days + OVERRUN);                 // number of sim days for data aggregation must be greater than OVERRUN + days discarded (from the beginning of sim)
    const size_t num_weeks = (par->runLength - discard_days - OVERRUN)/7;   // number of full weeks of data to be aggregated

    //vector<size_t> infected    = community->getNumNewlyInfected();
    vector<size_t> symptomatic = community->getNumNewlySymptomatic();
    vector<size_t> severe      = community->getNumNewlySevere();
    vector<size_t> dead        = community->getNumNewlyDead();

    // pair of num of primary infections starting this day, and mean num secondary infections they cause
    vector<pair<size_t, double>> R = community->getMeanNumSecondaryInfections();
    vector<size_t> Rt_incidence_tally(num_weeks, 0);

    vector<double> metrics(num_weeks*4, 0.0);
    for (size_t t = discard_days; t < discard_days + (7*num_weeks); ++t) {
        const size_t w = (t-discard_days)/7; // which reporting week are we in?
        metrics[w]                 += symptomatic[t];
        metrics[num_weeks + w]     += severe[t];
        metrics[2 * num_weeks + w] += dead[t];
        metrics[3 * num_weeks + w] += R[t].first > 0 ? R[t].first*R[t].second : 0;
        Rt_incidence_tally[w]      += R[t].first;
    }

    for (size_t w = 0; w < num_weeks; ++w) {
        metrics[3 * num_weeks + w] /= Rt_incidence_tally[w] > 0 ? Rt_incidence_tally[w] : 1.0;
    }
//metrics.resize(300);
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

//void calculate_reporting_ratios(Community* community) {
//    // this counts infections/cases/deaths that happen during the simulation,
//    // but not cases and deaths that are scheduled to happen after the last simulated day
//    const double cinf   = sum(community->getNumNewlyInfected());
//    const double ccase  = sum(community->getNumNewlySymptomatic());
//    const double csev   = sum(community->getNumNewlySevere());
//    const double ccrit  = sum(community->getNumNewlyCritical());
//    const double cdeath = sum(community->getNumNewlyDead());
//
//    cerr << "true infections, mild, severe, critical, deaths: " << cinf << " " << ccase << " " << csev  << " " << ccrit << " " << cdeath << endl;
//    cerr << "IFR, CFR: " << 100*cdeath/cinf << "%, " << 100*cdeath/ccase << "%" << endl;
//
//    // This is not a perfect way of calculating the case:death ratios, but it should be a reasonable approximation
//    const vector<vector<double>> RF = as_reported_fractions();
//    cerr << "\nReported fractions (asymp, mild, severe, crit, death [case:death ratio]):" << endl;
//    cerr << "wave 1: "; cerr_vector(RF[0]); cerr << " [" << (int) (RF[0][0]*cinf + RF[0][1]*ccase + RF[0][2]*csev + RF[0][3]*ccrit + RF[0][4]*cdeath)/(RF[0][4]*cdeath) << "] " << endl;
//    cerr << "wave 2: "; cerr_vector(RF[1]); cerr << " [" << (int) (RF[1][0]*cinf + RF[1][1]*ccase + RF[1][2]*csev + RF[1][3]*ccrit + RF[1][4]*cdeath)/(RF[1][4]*cdeath) << "] " << endl;
//    cerr << "wave 3: "; cerr_vector(RF[2]); cerr << " [" << (int) (RF[2][0]*cinf + RF[2][1]*ccase + RF[2][2]*csev + RF[2][3]*ccrit + RF[2][4]*cdeath)/(RF[2][4]*cdeath) << "] " << endl;
//
//    cerr << "\nIncidence by outcome:\n";
//    cerr << "\t ASYMPTOMATIC :\t" << community->getCumulIncidenceByOutcome(ASYMPTOMATIC) << endl;
//    cerr << "\t MILD :  \t" << community->getCumulIncidenceByOutcome(MILD) << endl;
//    cerr << "\t SEVERE :\t" << community->getCumulIncidenceByOutcome(SEVERE) << endl;
//    cerr << "\t CRITICAL :\t" << community->getCumulIncidenceByOutcome(CRITICAL) << endl;
//    cerr << "\t DEATH :\t" << community->getCumulIncidenceByOutcome(DEATH) << endl;
//}

vector<double> simulator(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par* mp = nullptr) {
    cerr << "rng seed: " << rng_seed << endl;
    gsl_rng_set(RNG, rng_seed);
    gsl_rng_set(VAX_RNG, rng_seed);
    gsl_rng_set(REPORTING_RNG, rng_seed);

//for (int i = 0; i < 1e6; ++i) { cerr << "RNG " << setprecision(40) << gsl_rng_uniform(RNG) << endl; } exit(10);
//for (int i = 0; i < 1e7; ++i) { cerr << "REPORTING_RNG " << gsl_rng_uniform(REPORTING_RNG) << endl; } exit(10);

    //gsl_rng_set(RNG, 1);
    //gsl_rng_set(VAX_RNG, 1);
    // initialize bookkeeping for run
    time_t start, end;
    time (&start);
    //const string process_id = report_process_id(args, serial, mp, start);
    //vector<double> abc_args(&args[0], &args[8]);
    vector<double> abc_args(args);
    //const size_t realization = 0; //(int) args[9];

    //const string process_id = report_process_id(abc_args, serial, start) + "." + to_string(realization);
    const string process_id = to_string(rng_seed);
    report_process_id(args, serial, GLOBAL_START_TIME, start);

    cerr << "SCENARIO " << rng_seed;
    for (auto _p: args) { cerr << " " << _p; } cerr << endl;

    // const size_t realization               = (size_t) args[0];
    const bool vaccine                     = 1; //(bool) args[1];
    const bool mutation                    = (bool) args[2];
    const size_t counterfactual_scenario   = 0; //(size_t) args[3];       // 0 = FL; 1 = VT; 2 = MS (should be set to 0 for all active strat work)
    const size_t dose_file                 = 0; //(size_t) args[4];       // 0 = state_based_counterfactual_doses.txt; 1 = active_vax_counterfactual_doses.txt; 2 =ring_vax_deployment_counterfactual_doses.txt
    const VacCampaignType active_vax_strat = (VacCampaignType) 0; //(size_t) args[5];       // 0 = none; 1 = ring vax; 2 = risk group vax; 3 = risk vax
    const bool quarantine_ctrl             = false; //(bool) args[6];     // 0 = off; 1 = on

    Parameters* par = define_simulator_parameters(args, rng_seed, serial, process_id);
    define_strain_parameters(par);

    vector<string> mutant_intro_dates = {};
    if (mutation) { mutant_intro_dates = {"2021-02-01", "2021-05-27", "2021-11-26"}; }

    Community* community = build_community(par);

    Vac_Campaign* vc = nullptr;
    community->setVac_Campaign(vc);

    par->immunityLeaky           = true;          // applies to both infection and vaccine immunity
    par->immunityWanes           = false;         // related to time-dep waning of protection, not waning of Ab levels
    par->seroPositivityThreshold = 0.0;

    // hanlde all vac campaign setup
    if (vaccine) {
        par->numVaccineDoses       = 3;             // total doses in series
        par->vaccineDoseInterval   = {21, 240};     // intervals between dose 1-->2 and dose 2-->3
        // par->vaccineTargetCoverage = 0.60;          // for healthcare workers only
        par->vaccine_dose_to_protection_lag = 10;   // number of days from vaccination to protection
        par->urgent_vax_dose_threshold = 1;         // the highest dose in series that will be administered in the active strategy

        const vector<string> vacFilenames = {"./state_based_counterfactual_doses.txt",
                                             "./active_vax_counterfactual_doses.txt",
                                             "./ring_vax_deployment_active_counterfactual_doses.txt",
                                             "./ring_vax_deployment_passive_counterfactual_doses.txt"};

        par->vaccinationFilename = vacFilenames[dose_file];

        // pooling will accumulate doses across age bins BUT NOT across doses
        bool pool_std_doses    = false;
        bool pool_urg_doses    = true;
        bool pool_all_doses    = false;

        // control whether to adjust to bin pops or total pop
        bool adjust_std_to_bin_pop = true;
        bool adjust_urg_to_bin_pop = false;

        vc = generateVac_Campaign(par, community, counterfactual_scenario, {pool_urg_doses, pool_std_doses, pool_all_doses}, adjust_std_to_bin_pop, adjust_urg_to_bin_pop);

        // parameter handling --- how do we want to handle setting these? I just set them here rather than use par
        vc->set_prioritize_first_doses(false);
        vc->set_flexible_queue_allocation(false);


        vc->set_reactive_vac_strategy(active_vax_strat);
        // vc->set_reactive_vac_dose_allocation(0.0);

        // the GROUPED_RISK_VACCINATION strategies starts alongside the GENERAL_CAMPAIGN
        // other active strategies start on May 1, 2021 (and if they require contact tracing, the check is called)
        int active_strat_start = active_vax_strat == GROUPED_RISK_VACCINATION
                                     ? vc->get_start_of_campaign(GENERAL_CAMPAIGN)
                                     : Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2021-05-01");
        if (vc->contact_tracing_required(active_vax_strat)) { assert(active_strat_start >= par->beginContactTracing); }
        vc->set_start_of_campaign(active_vax_strat, active_strat_start);

        vc->set_end_of_campaign(GENERAL_CAMPAIGN, par->runLength);
        vc->set_end_of_campaign(active_vax_strat, par->runLength);

        vector<int> min_ages(par->runLength, 5);
        vc->set_min_age(min_ages);       // needed for e.g. urgent vaccinations

        cerr << "Vaccination scenario:\n"
             << GENERAL_CAMPAIGN << "\n"
                                 << "\tduration " << vc->get_start_of_campaign(GENERAL_CAMPAIGN) << "--" << vc->get_end_of_campaign(GENERAL_CAMPAIGN) << "\n"
                                 << "\tnum doses " << par->numVaccineDoses << "\n"
                                 << "\tdose intervals "; cerr_vector(par->vaccineDoseInterval); cerr << "\n"
                                 << "\tdose protection lag " << par->vaccine_dose_to_protection_lag << "\n"
             << active_vax_strat   << "\n"
                                 << "\tduration " << vc->get_start_of_campaign(active_vax_strat) << "--" << vc->get_end_of_campaign(active_vax_strat) << "\n"
                                 << "\trequires contact tracing? " << boolalpha << vc->contact_tracing_required(active_vax_strat) << noboolalpha << "\n"
             << "other details"  << "\n"
                                 << "\tdose file " << par->vaccinationFilename << "\n"
                                 << "\tdose pooling (urg,std,all)? " << boolalpha << pool_urg_doses << ' ' << pool_std_doses << ' ' << pool_all_doses << noboolalpha << "\n"
                                 << "\tadj std doses to bin pop? " << boolalpha << adjust_std_to_bin_pop << noboolalpha << "\n"
                                 << "\tadj urg doses to bin pop? " << boolalpha << adjust_urg_to_bin_pop << noboolalpha << "\n"
                                 << "\tcontact tracing start " << par->beginContactTracing << "\n"
                                 << "\tself quarantining probs "; cerr_vector(par->quarantineProbability); cerr << "\n"
                                 << "\tself quarantining duration " << par->quarantineDuration << "\n" << endl;
    }
    // probability of self-quarantining for index cases and subsequent contacts
    if (quarantine_ctrl) {
        par->quarantineProbability = {0.9, 0.75, 0.5};
    } else {
        par->quarantineProbability = {0.0, 0.0, 0.0};
    }
    par->quarantineDuration = 10;

    // seed_epidemic(par, community, WILDTYPE);
    vector<string> plot_log_buffer = simulate_epidemic(par, community, process_id, mutant_intro_dates);//, social_contact_map);

    vector<double> cases(par->runLength, 0.0);
    vector<double> deaths(par->runLength, 0.0);

    for (auto person: community->getPeople()) {
        if (person->getNumNaturalInfections()) {
            vector<Infection*> infections = person->getInfectionHistory();
            for (Infection* inf: infections) {
                if (inf->getDetection()) {
                    const int inf_date = inf->getDetection()->reported_time;
                    if (inf_date < (int) cases.size()) {
                        cases[inf_date]++;
                        if (inf->fatal()) {
                            deaths[inf_date]++;
                        }
                    }
                }
            }
        }
    }

//    Date dummy_date(par);
//    for (unsigned int i = 0; i < cases.size(); ++i ) {
//        cerr << "cfr: " << i << " " << dummy_date.to_ymd() << " " << std::setprecision(4) << deaths[i] << " " << cases[i] << endl;
//        dummy_date.increment();
//    }

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
    string filename = "plot_log" + to_string(serial) + ".csv";
    write_daily_buffer(plot_log_buffer, process_id, filename, overwrite);
    stringstream ss;
    ss << "Rscript expanded_simvis.R " << serial;
    string cmd_str = ss.str();
    int retval = system(cmd_str.c_str());
    if (retval == -1) { cerr << "System call to `Rscript expanded_simvis.R` failed\n"; }
}

    time (&end);
    double dif = difftime (end,start);

    map<string, vector<vector<size_t>>> ages_by_outcome;
    ages_by_outcome["cases"].resize(2*par->runLength/7);
    ages_by_outcome["hosp"].resize(2*par->runLength/7);
    ages_by_outcome["deaths"].resize(2*par->runLength/7);

    vector<double> metrics = tally_counts(par, community, 0);
    //calculate_reporting_ratios(community);


    stringstream ss;
    ss << mp->mpi_rank << " end " << hex << process_id << " " << dec << dif << " ";

    for (auto i: args) ss << i << " ";
    for (auto i: metrics) ss << i << " ";
    ss << endl;

    string output = ss.str();
    fputs(output.c_str(), stderr);

    if (par->dump_simulation_data) {
        vector<string> tables = {
            "infection_history",
            "secondary_infections",
            "infection_detection",
            "vaccination_history",
            "age_bins",
            "doses_available",
            "doses_used"
        };
        generate_sim_data_db(par, community, serial, tables);
    }

    // if (vc)      { delete vc; }           // should this be here? MOVED INTO COMMUNITY DESTRUCTOR
    if (VAX_RNG) { gsl_rng_free(VAX_RNG); }
    delete community;
    delete par;

metrics = vector<double>(424);
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
//    if (not (argc == 3 or argc == 5 or argc == 6) ) {
//        usage();
//        exit(100);
//    }

    bool process_db = false;
    bool simulate_db = false;
    int buffer_size = 1;
    int requested_serial = -1;
    int requested_posterior_idx = -1;

    for (int i=2; i < argc;  i++ ) {
        if ( string_matches(argv[i], "--process") ) {
            process_db = true;
        } else if ( string_matches(argv[i], "--simulate") ) {
            simulate_db = true;
        } else if ( string_matches(argv[i], "-n" ) ) {
            buffer_size = atoi(argv[++i]);
        } else if ( string_matches(argv[i], "--serial" ) ) {
            requested_serial = atoi(argv[++i]);
        } else if ( string_matches(argv[i], "--posterior" ) ) {
            requested_posterior_idx = atoi(argv[++i]);
        } else if ( strcmp(argv[i], "--runLength" ) == 0 ) {
            TOTAL_DURATION = atoi(argv[++i]);
        } else if ( strcmp(argv[i], "--autotune" ) == 0 ) {
            autotune = true;
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

    delete abc; 
    return 0;
}
