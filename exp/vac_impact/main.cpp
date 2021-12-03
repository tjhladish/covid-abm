#include <chrono>
#include <unistd.h>
#include "AbcSmc.h"
#include "simulator.h"
#include <cstdlib>
#include "CCRC32.h"
#include "Utility.h"
#include <unordered_set>
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
const string vaccination_file = pop_dir + "/../fl_vac/fl_vac_v3.txt";

const int RESTART_BURNIN          = 0;
const int FORECAST_DURATION       = 700;
//const int FORECAST_DURATION       = 456;
const int OVERRUN                 = 14; // to get accurate Rt estimates near the end of the forecast duration
const bool RUN_FORECAST           = true;
const int TOTAL_DURATION          = RUN_FORECAST ? RESTART_BURNIN + FORECAST_DURATION + OVERRUN : RESTART_BURNIN;
//const size_t JULIAN_TALLY_DATE    = 146; // intervention julian date - 1
const size_t JULIAN_START_YEAR    = 2020;
//const double DEATH_UNDERREPORTING = 11807.0/20100.0; // FL Mar15-Sep5, https://www.nytimes.com/interactive/2020/05/05/us/coronavirus-death-toll-us.html

const int FL_POP                  = 21538187;   // as of 2020 census

vector<vector<double>> REPORTED_FRACTIONS;

const vector< vector<int> > VAX_AGE_BINS{ {12,14}, {15,24}, {25,34}, {35,44}, {45,54}, {55,64}, {65,74}, {75,84}, {85,NUM_AGE_CLASSES-1} };

const gsl_rng* VAX_RNG = gsl_rng_alloc(gsl_rng_mt19937);

//Parameters* define_simulator_parameters(vector<double> args, const unsigned long int rng_seed) {
Parameters* define_simulator_parameters(vector<double> /*args*/, const unsigned long int rng_seed, const unsigned long int serial, const string /*process_id*/) {
    Parameters* par = new Parameters();
    par->define_defaults();
    par->serial = serial;

    const float T = 0.0355;//0.022; // 0.0215 for the fl pop, 0.022 for the pseudo 1000k pop

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

    vector<double> seasonality;
    for (size_t day = 0; day < 366; ++day) {
        seasonality.push_back(1 - 0.15*sin((4*M_PI*((int) day-60))/366));
    }
    par->seasonality = seasonality;

    par->traceContacts = true; // needed for Rt calculation

    {
        // Create model for how outcome-dependent detection probabilities change over time
        //par->reportedFraction = {0.0, 0.01, 0.5, 0.8, 1.0};      // fraction of asymptomatic, mild, severe, critical, and deaths reported
        //par->reportedFraction = {0.0, 0.2, 0.75, 0.75, 0.75};      // fraction of asymptomatic, mild, severe, critical, and deaths reported
        //par->probFirstDetection = {0.0, 0.12, 0.55, 0.1, 0.01};      // probability of being detected while {asymp, mild, severe, crit, dead} if not detected previously

        const double RF_death_early = 0.8;//0.78; // overall probability of detecting death, at any point

        // probability of being detected while {asymp, mild, severe, crit, dead} if not detected previously
        vector<double> initial_vals = {0.0, 0.05, 0.7, 0.1};
        const double rho_death_ini  = 1.0 - (1.0 - RF_death_early)/((1.0 - initial_vals[0])*(1.0 - initial_vals[1])*(1.0 - initial_vals[2])*(1.0 - initial_vals[3]));
        initial_vals.push_back(rho_death_ini);

        vector<double> mid_vals   = {0.02, 0.5, 0.1, 0.1};
        const double rho_death_mid  = 1.0 - (1.0 - RF_death_early)/((1.0 - mid_vals[0])*(1.0 - mid_vals[1])*(1.0 - mid_vals[2])*(1.0 - mid_vals[3]));
        mid_vals.push_back(rho_death_mid);

        const double RF_death_late = 1.0;//0.78; // overall probability of detecting death, at any point
        vector<double> final_vals   = {0.1, 0.9, 0.9, 0.0};
        const double rho_death_fin  = 1.0 - (1.0 - RF_death_late)/((1.0 - final_vals[0])*(1.0 - final_vals[1])*(1.0 - final_vals[2])*(1.0 - final_vals[3]));
        final_vals.push_back(rho_death_fin);

        cerr << "death init, mid, fin: " << rho_death_ini << " " << rho_death_mid << " " << rho_death_fin << endl;

        const int isd1 = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2020-06-01");
        const int isd2 = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2020-10-01");

        const vector<int> inflection1_sim_day = {isd1, isd1, isd1, isd1, isd1};
        const vector<int> inflection2_sim_day = {isd2, isd2, isd2, isd2, isd2};
        const vector<double> slopes = {0.1, 0.1, 0.1, 0.1, 0.1}; // sign is determined based on initial/final values

        par->createDetectionModel(initial_vals, mid_vals, final_vals, inflection1_sim_day, inflection2_sim_day, slopes, slopes);

        vector<double> reported_frac_init  = par->toReportedFraction(initial_vals);
        vector<double> reported_frac_mid   = par->toReportedFraction(mid_vals);
        vector<double> reported_frac_final = par->toReportedFraction(final_vals);
        cerr_vector(reported_frac_init); cerr << endl;  REPORTED_FRACTIONS.push_back(reported_frac_init);
        cerr_vector(reported_frac_mid); cerr << endl;   REPORTED_FRACTIONS.push_back(reported_frac_mid);
        cerr_vector(reported_frac_final); cerr << endl; REPORTED_FRACTIONS.push_back(reported_frac_final);

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

    vector<TimeSeriesAnchorPoint> ap = { // tuning for whole FL model
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
    par->pathogenicityReduction = 0.0;                      // to be fit; fraction of infections missed in pathogenicity studies
                                                            // used for interpreting input data, *not within simulation*
//    par->susceptibilityCorrection = 1.0;
    par->define_susceptibility_and_pathogenicity();

    //par->daysImmune = 730; // changing this to be a function call
    par->VES = {{WILDTYPE, {0.0}}};

    //par->hospitalizedFraction = 0.25; // fraction of cases assumed to be hospitalized
    par->probInitialExposure = {4.0e-04};
    //par->probDailyExposure   = vector(120, 2.0e-05);        // introductions are initially lower
    //par->probDailyExposure.resize(par->runLength, 2.0e-04); // and then pick up after ~ 4 months
    par->probDailyExposure   = {1.0e-04};        // introductions are initially lower

    par->populationFilename       = pop_dir    + "/population-"         + SIM_POP + ".txt";
    par->comorbidityFilename      = pop_dir    + "/comorbidity-"        + SIM_POP + ".txt";
    par->locationFilename         = pop_dir    + "/locations-"          + SIM_POP + ".txt";
    par->networkFilename          = pop_dir    + "/network-"            + SIM_POP + ".txt";
    par->publicActivityFilename   = pop_dir    + "/public-activity-"    + SIM_POP + ".txt";

    return par;
}

void parseVaccineFile(string vaccinationFilename, const Parameters* par, Community* community, Vac_Campaign* vc, set<Person*>& scheduled_people, map<int, int>& sch_hcw_by_age) {
    // ratio of synthpop to FL pop
    const double pop_ratio = (double)community->getNumPeople()/FL_POP;

    // temporary data structure to hold calculated doses used before setting up the vac_campaign
    vector< vector<int> > doses_available;
    doses_available.resize(par->runLength, vector<int>(NUM_OF_VACCINE_ALLOCATION_TYPES));

    // create map of vectors of unvaccinated people by age bin (starts with all people) with key = min bin age
    map<int, vector<Person*> > unscheduled_people;
    // for(auto bin : VAX_AGE_BINS){ unscheduled_people[bin[0]] = vector<Person*>(); }

    // for each age bin...
    for(auto v : VAX_AGE_BINS) {
        const int age_min = v[0], age_max = v[1];
        // and for each age in the bin...
        for(int age = age_min; age <= age_max; ++age) {
            // grab all people of that age ONLY if they are not a HCW
            for(Person* p : community->getAgeCohort(age)) { if(not p->isHCW()) { unscheduled_people[age_min].push_back(p); } }
        }
        // shuffle age bin after all people are added
        gsl_ran_shuffle(VAX_RNG, unscheduled_people[age_min].data(), unscheduled_people[age_min].size(), sizeof(Person*));
    }

    // check that vaccinationFilename exists and can be opened
    ifstream iss(vaccinationFilename);
    string buffer;
    istringstream line;

    if (!iss) {
        cerr << "ERROR: vaccination file " << vaccinationFilename << " not found." << endl;
        exit(-1);
    }

    // variables for vaccination file processing
    string age_grp, end_of_week_date, prev_wk = "0000-00-00";
    size_t mrna_dose_1, mrna_dose_2, jj_doses;
    size_t bin_min, bin_max;

    // saves the final vaccination rates for each age bin for use in projected vaccination
    map< int, vector<double> > last_known_vac_rate;
    size_t last_known_revac_doses = 0;

    size_t last_day_of_data = 0;
    while (getline(iss, buffer)) {
        line.clear();
        line.str(buffer);

        // vector to hold all scheduled people (regardless of age) for a given week (for shuffling purposes)
        vector<Person*> people_to_be_scheduled;

        if (line >> age_grp >> end_of_week_date >> mrna_dose_1 >> mrna_dose_2 >> jj_doses >> bin_min >> bin_max) {
            // calculate and save vaccination rates into last_known_vac_rate
            const double rate_denom = unscheduled_people[bin_min].size() * 7.0;
            assert(rate_denom > 0);
            last_known_vac_rate[bin_min] = {mrna_dose_1/rate_denom, mrna_dose_2/rate_denom, jj_doses/rate_denom };    // used to project future vaccination rates

            // keeps track of which week is being processed
            // if the current week has changed, shuffle the previous week's people and schedule
            if(end_of_week_date != prev_wk) {
                if(people_to_be_scheduled.size() != 0){ gsl_ran_shuffle(VAX_RNG, people_to_be_scheduled.data(), people_to_be_scheduled.size(), sizeof(Person*)); }
                for(Person* p : people_to_be_scheduled) { vc->schedule_vaccination(p); }    // revaccinations automatically handled in Community::vaccinate()
                people_to_be_scheduled.clear();
                prev_wk = end_of_week_date;
            }

            // tally daily dose totals for this age bin
            const double daily_total_first_doses = (mrna_dose_1 + jj_doses)/7.0;
            // const double daily_total_doses = daily_total_first_doses + (mrna_dose_2/7.0);

            // binomial distribution parameters
            const double sch_binom_np = (daily_total_first_doses) * pop_ratio;        // np = total number of doses administered to this age group adjusted for synthpop size

            // use date (which is the final day of an week of reported data) from file to cycle through the days in the matching simulation week
            const size_t end_of_week = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, end_of_week_date);
            for(size_t day = end_of_week-6; day <= end_of_week; ++day) {
                const size_t revacDay = day+par->vaccineDoseInterval;
                last_day_of_data = max(day, last_day_of_data);
                if(day >= par->runLength) { continue; }
                // binomial distribution parameters
                const size_t sch_binom_n = unscheduled_people[bin_min].size();      // n = number of unscheduled people remaining in age bin
                const double sch_binom_p = sch_binom_np/sch_binom_n;                // p = probability of selecting someone from the age group to be vaccinated

                // select number of people to vaccinate for this week from this age group using binomial distribution
                const size_t binom_sample = gsl_ran_binomial(VAX_RNG, sch_binom_p, sch_binom_n);
                const size_t num_to_schedule = max(0, (int) binom_sample - sch_hcw_by_age[bin_min]);
                sch_hcw_by_age[bin_min]      = max(0, sch_hcw_by_age[bin_min] - (int) binom_sample);

                // schedule proper number of vaccinations from proper age group and remove from unscheduled data structure
                for(size_t num_scheduled = 0; num_scheduled < num_to_schedule; ++num_scheduled){
                    Person* p = unscheduled_people[bin_min].back();
                    // set.insert.second returns bool (true if inserted, false if not) --- this keeps track of who has been scheduled and prevents double-scheduling
                    if((scheduled_people.insert(p).second)){ people_to_be_scheduled.push_back(p); }
                    unscheduled_people[bin_min].pop_back();
                }

                // // tally doses used for that day for urgent allocation (reactive strategy) adjusted for synthpop size
                // const size_t emp_daily_urgent_doses = vc->get_reactive_vac_strategy() != NUM_OF_REACTIVE_VAC_STRATEGY_TYPES ?
                //                           (size_t) round(vc->get_reactive_vac_dose_allocation() * daily_total_doses * pop_ratio) :
                //                           0;

                // // tally doses used for that week for standard allocation (general strategy) adjusted for synthpop size
                // const size_t emp_daily_standard_doses = vc->get_reactive_vac_strategy() != NUM_OF_REACTIVE_VAC_STRATEGY_TYPES ?
                //                             (size_t) round((1 - vc->get_reactive_vac_dose_allocation()) * daily_total_doses * pop_ratio) :
                //                             (size_t) round(daily_total_doses * pop_ratio);

                const size_t emp_daily_standard_first_doses = (size_t) round(daily_total_first_doses * pop_ratio);
                // doses_available.at(day)[URGENT_ALLOCATION]   += emp_daily_urgent_doses;
                doses_available.at(day)[STANDARD_ALLOCATION] += emp_daily_standard_first_doses;
                if(revacDay < par->runLength) {
                    doses_available.at(revacDay)[STANDARD_ALLOCATION] += emp_daily_standard_first_doses;
                    last_known_revac_doses = doses_available.at(revacDay)[STANDARD_ALLOCATION];
                }
            }
        }
        // schedule any remaining people
        if(people_to_be_scheduled.size() != 0){ gsl_ran_shuffle(VAX_RNG, people_to_be_scheduled.data(), people_to_be_scheduled.size(), sizeof(Person*)); }
        for(Person* p : people_to_be_scheduled) { vc->schedule_vaccination(p); }
        people_to_be_scheduled.clear();
    }
    iss.close();

    // PROJECTED VACCINATION TO 2021-08-31
    // at same rate as 2021-05-29, vaccinate people until 2021-08-31
    if (par->runLength > last_day_of_data + 1) {
        vector<Person*> people_to_be_scheduled;

        // calculate the total number of unscheduled people left on 2021-05-29
        size_t num_current_unsch = 0;
        for(auto const& [bin_min, unsch_grp] : unscheduled_people) { num_current_unsch += unsch_grp.size(); }

        for(size_t day = last_day_of_data+1; day < par->runLength; ++day) {
            const size_t revacDay = day+par->vaccineDoseInterval;
            // const size_t proj_daily_urgent_doses =  * ;
            // const size_t proj_daily_standard_doses = doses_available.at(may29)[STANDARD_ALLOCATION];
            // adjust available doses by the ratio of unscheduled people on day to unscheduled people on 2021-05-29
            // results in diminishing number of doses available without causing dose shortages
            double proj_dose_adj = 0;
            for(auto bin : VAX_AGE_BINS) {
                const size_t bin_min = bin[0];
                const double sch_binom_p = pop_ratio*(last_known_vac_rate[bin_min][0] + last_known_vac_rate[bin_min][2]);

                // binomial distribution parameters
                const size_t sch_binom_n = unscheduled_people[bin_min].size();      // n = number of unscheduled people remaining in age bin
                proj_dose_adj += sch_binom_n;
                // select number of people to vaccinate for this week from this age group using binomial distribution
                // POSSIBLE TODO:  add binom_sample, sch_hcw_by_age adjustment like how it is above
                const size_t num_to_schedule = gsl_ran_binomial(VAX_RNG, sch_binom_p, sch_binom_n);
                // schedule proper number of vaccinations from proper age group and remove from unscheduled data structure
                for(size_t num_scheduled = 0; num_scheduled < num_to_schedule; ++num_scheduled){
                    Person* p = unscheduled_people[bin_min].back();
                    // set.insert.second returns bool (true if inserted, false if not) --- this keeps track of who has been scheduled and prevents double-scheduling
                    if((scheduled_people.insert(p).second)){ people_to_be_scheduled.push_back(p); }
                    unscheduled_people[bin_min].pop_back();
                }
            }
            if(people_to_be_scheduled.size() != 0){ gsl_ran_shuffle(VAX_RNG, people_to_be_scheduled.data(), people_to_be_scheduled.size(), sizeof(Person*)); }
            for(Person* p : people_to_be_scheduled) { vc->schedule_vaccination(p); }
            people_to_be_scheduled.clear();

            // vaccinate same fraction of unvacinated people each day
            // doses_available.at(day)[URGENT_ALLOCATION]   += round(proj_dose_adj*doses_available.at(last_day_of_data)[URGENT_ALLOCATION] / num_current_unsch);
            doses_available.at(day)[STANDARD_ALLOCATION] += round(proj_dose_adj*last_known_revac_doses / num_current_unsch);
            if(revacDay < par->runLength) { doses_available.at(revacDay)[STANDARD_ALLOCATION] += round(proj_dose_adj*last_known_revac_doses / num_current_unsch); }
        }
    }
    vc->set_doses_available(doses_available);

}

Vac_Campaign* generateVac_Campaign(string vaccinationFilename, const Parameters* par, Community* community) {
    // keeps track of who has been scheduled to prevent scheduling the same person twice
    set<Person*> scheduled_people;

    // create a new Vac_Campaign
    Vac_Campaign* vc = new Vac_Campaign();

    map<int, int> sch_hcw_by_age;
    map<int, int> age_to_bin_map;
    for(auto& bin : VAX_AGE_BINS) {
        for(int i = bin[0]; i <= bin[1]; ++i) { age_to_bin_map[i] = bin[0]; }
    }

    // add all healthcare workers (hospital + nursing home workers) to standard queue BEFORE all other people are scheduled
    for(Person* p : community->getPeople()) {
        // get work location ID + check if workplace is a hospital or nursing home
        // if so, schedule vaccination and add to scheduled_people
        if(p->isHCW()) {
            // set.insert.second returns bool (true if inserted, false if not) --- this keeps track of who has been scheduled and prevents double-scheduling
            if((gsl_rng_uniform(VAX_RNG) < par->vaccineTargetCoverage) and (scheduled_people.insert(p).second)){
                assert(p->getAge() >= 12);
                vc->schedule_vaccination(p);
                sch_hcw_by_age[age_to_bin_map.at(p->getAge())]++;
            }
        }
    }

    parseVaccineFile(vaccinationFilename, par, community, vc, scheduled_people, sch_hcw_by_age);
    community->setVac_Campaign(vc);

    return vc;
}

// Take a list of values, return original indices sorted by value
vector<int> ordered(vector<int> const& values) {

    vector<pair<int,int> > pairs(values.size());
    for(size_t pos=0; pos<values.size(); pos++) {
        pairs[pos] = make_pair(values[pos],pos);
    }

    //bool comparator ( const mypair& l, const mypair& r) { return l.first < r.first; }
    std::sort( pairs.rbegin(), pairs.rend() ); // sort greatest to least
    vector<int> indices(values.size());
    for(size_t i=0; i < pairs.size(); i++) indices[i] = pairs[i].second;

    return indices;
}


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
    assert((int) par->runLength >= discard_days + OVERRUN);                 // number of sim days for data aggregation must be greater than OVERRUN + days discarded (from the beginning of sim)
    const size_t num_weeks = (par->runLength - discard_days - OVERRUN)/7;   // number of full weeks of data to be aggregated

    //vector<size_t> infected    = community->getNumNewlyInfected();
    //vector< vector<int> > severe      = community->getNumSevereCases();
    vector<size_t> symptomatic = community->getNumNewlySymptomatic();
    vector<size_t> dead        = community->getNumNewlyDead();

    // pair of num of primary infections starting this day, and mean num secondary infections they cause
    vector<pair<size_t, double>> R = community->getMeanNumSecondaryInfections();
    vector<size_t> Rt_incidence_tally(num_weeks, 0);

    vector<double> metrics(num_weeks*3, 0.0);
    for (size_t t = discard_days; t < discard_days + (7*num_weeks); ++t) {
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

void calculate_reporting_ratios(Community* community) {
    // this counts infections/cases/deaths that happen during the simulation,
    // but not cases and deaths that are scheduled to happen after the last simulated day
    const double cinf   = sum(community->getNumNewlyInfected());
    const double ccase  = sum(community->getNumNewlySymptomatic());
    const double csev   = sum(community->getNumNewlySevere());
    const double ccrit  = sum(community->getNumNewlyCritical());
    const double cdeath = sum(community->getNumNewlyDead());

    cerr << "true infections, mild, severe, critical, deaths: " << cinf << " " << ccase << " " << csev  << " " << ccrit << " " << cdeath << endl;
    cerr << "IFR, CFR: " << 100*cdeath/cinf << "%, " << 100*cdeath/ccase << "%" << endl;

    // This is not a perfect way of calculating the case:death ratios, but it should be a reasonable approximation
    const vector<vector<double>> RF = REPORTED_FRACTIONS;
    cerr << "\nReported fractions (asymp, mild, severe, crit, death [case:death ratio]):" << endl;
    cerr << "wave 1: "; cerr_vector(RF[0]); cerr << " [" << (int) (RF[0][0]*cinf + RF[0][1]*ccase + RF[0][2]*csev + RF[0][3]*ccrit + RF[0][4]*cdeath)/(RF[0][4]*cdeath) << "] " << endl;
    cerr << "wave 2: "; cerr_vector(RF[1]); cerr << " [" << (int) (RF[1][0]*cinf + RF[1][1]*ccase + RF[1][2]*csev + RF[1][3]*ccrit + RF[1][4]*cdeath)/(RF[1][4]*cdeath) << "] " << endl;
    cerr << "wave 3: "; cerr_vector(RF[2]); cerr << " [" << (int) (RF[2][0]*cinf + RF[2][1]*ccase + RF[2][2]*csev + RF[2][3]*ccrit + RF[2][4]*cdeath)/(RF[2][4]*cdeath) << "] " << endl;

    cerr << "\nIncidence by outcome:\n";
    cerr << "\t ASYMPTOMATIC :\t" << Community::_cumulIncByOutcome[ASYMPTOMATIC] << endl;
    cerr << "\t MILD :  \t" << Community::_cumulIncByOutcome[MILD] << endl;
    cerr << "\t SEVERE :\t" << Community::_cumulIncByOutcome[SEVERE] << endl;
    cerr << "\t CRITICAL :\t" << Community::_cumulIncByOutcome[CRITICAL] << endl;
    cerr << "\t DEATH :\t" << Community::_cumulIncByOutcome[DEATH] << endl;
}

vector<double> simulator(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par* mp = nullptr) {
    cerr << "rng seed: " << rng_seed << endl;
    gsl_rng_set(RNG, rng_seed); // seed the rng using sys time and the process id
    gsl_rng_set(VAX_RNG, rng_seed); // seed the rng using sys time and the process id
    // initialize bookkeeping for run
    time_t start, end;
    time (&start);
    //const string process_id = report_process_id(args, serial, mp, start);
    //vector<double> abc_args(&args[0], &args[8]);
    vector<double> abc_args(args);
    //const size_t realization = 0; //(int) args[9];

    //const string process_id = report_process_id(abc_args, serial, start) + "." + to_string(realization);
    const string process_id = to_string(rng_seed);
    report_process_id(args, serial, start);

    cerr << "SCENARIO " << rng_seed;
    for (auto _p: args) { cerr << " " << _p; } cerr << endl;

    Parameters* par = define_simulator_parameters(args, rng_seed, serial, process_id);
    const bool vaccine             = (bool) args[0];
    // const size_t realization    = (size_t) args[1];
    const bool mutation          = (bool) args[2];
    vector<string> mutant_intro_dates = {};
    if (mutation) { mutant_intro_dates = {"2021-02-01", "2021-05-27"}; };   // extra semicolon?

    Community* community = build_community(par);
    Vac_Campaign* vc = nullptr;
    community->setVac_Campaign(vc);

    if (vaccine) {
    // TODO -- Schedule health care workers first, then either use actual age dist data, or sample from over 65, over 16, over 12, etc.
        //double target_coverage  = coverage;

        // below was previous vaccination system
        /*
        double catchup_coverage = 0.6;
        const int target = 16;
        const int senior_threshold = 65;
        size_t senior_vac_campaign_duration = vac_rate == SLOW ? 86 : 35; // 0.12% per day vs 0.3% per day for [65, max_age]
        size_t general_vac_campaign_duration = vac_rate == SLOW ? 257 : 103; // same, for those [16, 64]
        const size_t senior_vac_sim_day  = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, "2021-01-08"); // jan 1 + 7 days for efficacy
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
        */
        const double x_alpha_1 = 1;//0.529; // calculated using quadratic formula: (VES*VEP)x^2 - (VES+VEP)x + VESP_alpha = 0, where VES and VEP are for wildtype
        const double x_alpha_2 = 1;//0.948;

        const double x_delta_1 = 1;//0.338;
        const double x_delta_2 = 1;//0.843;

        par->VES                   = {{WILDTYPE, {0.40, 0.80}}, {B_1_1_7, {0.40*x_alpha_1, 0.80*x_alpha_2}}, {B_1_617_2, {0.40*x_delta_1, 0.80*x_delta_2}}}; // reduce for delta
        //par->VES                   = {{WILDTYPE, {0.40, 0.80}}, {B_1_1_7, {0.40*x_alpha_1, 0.80*x_alpha_2}}, {B_1_617_2, {0.05, 0.52}}}; // reduce for delta
        par->VES_NAIVE             = par->VES;
        par->VEP                   = {{WILDTYPE, {0.67, 0.75}}, {B_1_1_7, {0.67*x_alpha_1, 0.75*x_alpha_2}}, {B_1_617_2, {0.67*x_delta_1, 0.75*x_delta_2}}}; // reduce for delta
        //par->VEP                   = {{WILDTYPE, {0.67, 0.75}}, {B_1_1_7, {0.67*x_alpha_1, 0.75*x_alpha_2}}, {B_1_617_2, {0.29, 0.75}}}; // reduce for delta
        par->VEH                   = {{WILDTYPE, {0.9, 1.0}},   {B_1_1_7, {0.9, 1.0}}, {B_1_617_2, {0.9, 1.0}}};
        par->VEI                   = {{WILDTYPE, {0.4, 0.8}},   {B_1_1_7, {0.2, 0.5}}, {B_1_617_2, {0.1, 0.1}}};
        par->VEF                   = {{WILDTYPE, {0.0, 0.0}},   {B_1_1_7, {0.0, 0.0}}, {B_1_617_2, {0.0, 0.0}}}; // effect likely captured by VEH
        par->vaccineLeaky          = true;
        par->numVaccineDoses       = 2;
        par->vaccineDoseInterval   = 21;
        par->vaccineTargetCoverage = 0.60;  // for healthcare workers only
        par->vaccine_dose_to_protection_lag = 10;

        vc = generateVac_Campaign(vaccination_file, par, community);

        // parameter handling --- how do we want to handle setting these? I just set them here rather than use par
        vc->set_prioritize_first_doses(false);
        vc->set_flexible_queue_allocation(false);
        vc->set_reactive_vac_strategy(NUM_OF_VAC_CAMPAIGN_TYPES);
        vc->set_reactive_vac_dose_allocation(0.0);
        vector<int> min_ages(par->runLength, 12);
        vc->set_min_age(min_ages);       // needed for e.g. urgent vaccinations
    }

    seed_epidemic(par, community, WILDTYPE);
    vector<string> plot_log_buffer = simulate_epidemic(par, community, process_id, mutant_intro_dates);

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

    vector<double> metrics = tally_counts(par, community, 0);
    //calculate_reporting_ratios(community);

    stringstream ss;
    ss << mp->mpi_rank << " end " << hex << process_id << " " << dec << dif << " ";

    for (auto i: args) ss << i << " ";
    for (auto i: metrics) ss << i << " ";
    ss << endl;

    string output = ss.str();
    fputs(output.c_str(), stderr);

    if(vc) { delete vc; }           // should this be here?
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
//        } else if ( strcmp(argv[i], "--runLength" ) == 0 ) {
//            TOTAL_DURATION = atoi(argv[++i]);
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
