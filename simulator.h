#include <cstdlib>
#include <cstring>
#include <climits>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Parameters.h"
#include "Person.h"
#include "Location.h"
#include "Community.h"
#include "Utility.h"
#include "sys/stat.h"
#include "Vac_Campaign.h"
#include <utility>

using namespace covid::standard;
using namespace covid::util;

enum IncidenceReportingType {
    INTRO_INF,
    TOTAL_INF,
    INTRO_CASE,
    TOTAL_CASE,
    INTRO_SEVERE,
    TOTAL_SEVERE,
    INTRO_CRITICAL,
    TOTAL_CRITICAL,
    INTRO_DEATH,
    TOTAL_DEATH,
    NUM_OF_INCIDENCE_REPORTING_TYPES
};

enum PrevalenceReportingType {
    INTRO_INF_PREV,
    TOTAL_INF_PREV,
    INTRO_CASE_PREV,
    TOTAL_CASE_PREV,
    INTRO_SEVERE_PREV,
    TOTAL_SEVERE_PREV,
    INTRO_CRITICAL_PREV,
    TOTAL_CRITICAL_PREV,
    INTRO_DEATH_PREV,
    TOTAL_DEATH_PREV,
    NUM_OF_PREVALENCE_REPORTING_TYPES
};


gsl_rng* RNG = gsl_rng_alloc(gsl_rng_taus2);
gsl_rng* REPORTING_RNG = gsl_rng_alloc(gsl_rng_mt19937);

// Predeclare local functions
Community* build_community(const Parameters* par);
void seed_epidemic(const Parameters* par, Community* community, StrainType strain);
void write_immunity_file(const Community* community, const string label, string filename, int runLength);
void write_immunity_by_age_file(const Community* community, const int year, string filename="");
void write_daily_buffer( vector<string>& buffer, const string process_id, string filename="", bool overwrite=false);

Community* build_community(const Parameters* par) {
    Date* date = new Date(par);
    Community* community = new Community(par, date);
    Person::setPar(par);
cerr << "Reading locations ... ";
    if (!community->loadLocations(par->locationFilename, par->networkFilename)) {
        cerr << "ERROR: Could not load locations" << endl;
        exit(-1);
    }
cerr << "done.\n";
cerr << "Reading population ... ";
    if (!community->loadPopulation(par->populationFilename, par->comorbidityFilename, par->publicActivityFilename)) {
        cerr << "ERROR: Could not load population" << endl;
        exit(-1);
    }

cerr << "done.\n"; //  Now sleeping for 20s so ram usage can be checked.\n";
//sleep(20);
    if (par->abcVerbose) {
        cerr << community->getNumPeople() << " people" << endl;
    }

    /*if (!par->bSecondaryTransmission) {
        community->setNoSecondaryTransmission();
    }*/

    return community;
}

// DEPRECATED
// Community* deep_copy_community(const Parameters* par) {
//     Date* date = new Date(par);
//     Community* community = new Community(par, date);
//     Person::setPar(par);
//     // deep copy locations
//     // deep copy population
//     return community;
// }


void seed_epidemic(const Parameters* par, Community* community, StrainType strain) {
    // epidemic may be seeded with initial exposure OR initial infection
    bool attempt_initial_infection = true;
    // Normal usage, to simulate epidemic
    const size_t pop_size = community->getNumPeople();
    if (par->numInitialExposed > 0) {
        attempt_initial_infection = false;
        for (size_t i=0; i<par->numInitialExposed; i++) {
            community->infect(gsl_rng_uniform_int(RNG, pop_size), strain);
        }
    } else if (par->probInitialExposure > 0.0) {
        // determine how many people are exposed
        size_t k = gsl_ran_binomial(RNG, par->probInitialExposure, pop_size);
        vector<size_t> pids(community->getNumPeople());
        iota(pids.begin(), pids.end(), 0);
        vector<size_t> exposed_group = choose_k(RNG, pids, k);

        size_t inf_ct = 0;
        for (auto pid: exposed_group) {
            Infection* inf = community->infect(pid, strain);
            if (inf) { inf_ct++; }
        }
        cerr << "pop size, sampled size, infected size: " << pop_size << ", " << k << ", " << inf_ct << endl;
    }

    if (attempt_initial_infection) {
        // Useful for estimating R0
        if(par->numInitialInfected > 0) {
            int count = community->getNumInfected(0);

            // must infect initialInfected persons -- this bit is mysterious
            while (community->getNumInfected(0) < count + par->numInitialInfected) {
                community->infect(gsl_rng_uniform_int(RNG, pop_size), strain);
            }
        }
    }
    return;
}


void _aggregator(map<string, vector<int> >& periodic_incidence, string key) {
    for (unsigned int i = 0; i < periodic_incidence["daily"].size(); ++i) periodic_incidence[key][i] += periodic_incidence["daily"][i];
}


void _reporter(stringstream& ss, map<string, vector<int> > &periodic_incidence, vector<int> &periodic_prevalence, const Parameters* par, const string process_id, const string label, const int value, string key) {
        ss << process_id << dec << " " << par->serial << label << value << " ";
        for (auto v: periodic_incidence[key]) ss << v << " ";
        if(key=="daily") for (auto v: periodic_prevalence) ss << v << " ";
}


void periodic_output(const Parameters* par, map<string, vector<int> > &periodic_incidence, vector<int> &periodic_prevalence, const Date* date, const string process_id, vector<int>& epi_sizes) {
    stringstream ss;

    if (par->dailyOutput) {
        _reporter(ss, periodic_incidence, periodic_prevalence, par, process_id, " day: ", date->day(), "daily"); ss << endl;
    }

    vector<int> dummy;
     if (par->periodicOutput) {
        _aggregator(periodic_incidence, "n_day");
        const int n = par->periodicOutputInterval;
        if (date->endOfPeriod(n)) {
            _reporter(ss, periodic_incidence, dummy, par, process_id, " " + to_string(n) + "_day: ", date->nDayPeriod(n), "n_day"); ss << endl;
            periodic_incidence["n_day"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
        }
    }

    if (par->weeklyOutput) {
        _aggregator(periodic_incidence, "weekly");
        if (date->endOfWeek()) {
            _reporter(ss, periodic_incidence, dummy, par, process_id, " week: ", date->week(), "weekly"); ss << endl;
            periodic_incidence["weekly"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
        }
    }

    if (par->monthlyOutput) {
        _aggregator(periodic_incidence, "monthly");
        if (date->endOfMonth()) {
            _reporter(ss, periodic_incidence, dummy, par, process_id, " month: ", date->julianMonth(), "monthly"); ss << endl;
            periodic_incidence["monthly"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
        }
    }

    // handle several things that happen yearly
    _aggregator(periodic_incidence, "yearly");
    if (date->endOfYear()) {
        if (par->abcVerbose) {
            cout << process_id << dec << " " << par->serial << " T: " << date->day() << " annual: ";
            for (auto v: periodic_incidence["yearly"]) { cout << v << " "; } cout << endl;
        }

        epi_sizes.push_back(periodic_incidence["yearly"][2]);

        if (par->yearlyOutput) { _reporter(ss, periodic_incidence, dummy, par, process_id, " year: ", date->year(), "yearly"); ss << endl; }
        periodic_incidence["yearly"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
    }

    periodic_incidence["daily"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
    periodic_prevalence = vector<int>(NUM_OF_PREVALENCE_REPORTING_TYPES, 0);
    string output = ss.str();
    //fputs(output.c_str(), stderr);
    fputs(output.c_str(), stdout);
}


int seed_epidemic(const Parameters* par, Community* community, const Date* date, vector<double> strain_weights) {
    int introduced_infection_ct = 0;
    const int numperson = community->getNumPeople();
    const size_t dailyExposedIdx = date->day() % par->probDailyExposure.size();
    const double intro_rate_multiplier = *date > "2021-06-15" ? 2.0 : 1.0;
    const double expected_num_exposed = intro_rate_multiplier * par->probDailyExposure[dailyExposedIdx] * numperson;
    if (expected_num_exposed > 0) {
        assert(expected_num_exposed <= numperson);
        const int num_exposed = gsl_ran_poisson(RNG, expected_num_exposed);
        for (int i=0; i<num_exposed; i++) {
            // gsl_rng_uniform_int returns on [0, numperson-1]
            int transmit_to_id = gsl_rng_uniform_int(RNG, numperson);
            if (community->infect(transmit_to_id, (StrainType) weighted_choice(RNG, strain_weights))) {
                introduced_infection_ct++;
            }
        }
    }
    return introduced_infection_ct;
}


map<string, vector<int> > construct_tally() {
    // { introductions, local transmission, total, case, severe}
    map<string, vector<int> > periodic_incidence { {"daily", vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES,0)},
                                                   {"n_day", vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES,0)},
                                                   {"weekly", vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES,0)},
                                                   {"monthly", vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES,0)},
                                                   {"yearly", vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES,0)} };
    return periodic_incidence;
}


template<class T>
size_t tally_decreases(const vector<T> &vals) {
    size_t hits = 0;
    for (size_t i = 1; i < vals.size(); ++i) {
        hits += vals[i-1] < vals[i];
    }
    return hits;
}


/*
Date,rcase,rdeath,death_incd
2020-03-02,2,0,0
2020-03-03,2,0,0
2020-03-04,4,0,0
*/
// parse empirical data file and store in a map with key of sim_day and value of a vector of cases and deaths
map<size_t, vector<int>> parse_emp_data_file(const Parameters* par, const string emp_data_file) {
    vector< vector<string> > emp_data = read_2D_vector_file(emp_data_file, ',');
    map<size_t, vector<int>> recast_emp_data;

    for (vector<string> &v : emp_data) {
        if (v[0] == "Date") { continue; }
        int rcase = stoi(v[1]), rdeath = stoi(v[3]);
        recast_emp_data[Date::to_sim_day(par->startJulianYear, par->startDayOfYear, v[0])] = {rcase, rdeath};
    }
    return recast_emp_data;
}

// encapsulates the data structures that are important to cache in order to minimize code footprint
class SimulationLedger {
public:
    SimulationLedger() = default;
    SimulationLedger(const SimulationLedger &o) {
        epi_sizes           = o.epi_sizes;
        daily_output_buffer = o.daily_output_buffer;
        periodic_incidence  = o.periodic_incidence;
        periodic_prevalence = o.periodic_prevalence;
        plot_log_buffer     = o.plot_log_buffer;
        strains             = o.strains;
    };
    virtual ~SimulationLedger() {};

    vector<int> epi_sizes;
    vector<string> daily_output_buffer;
    map<string, vector<int> > periodic_incidence;
    vector<int> periodic_prevalence;
    vector<string> plot_log_buffer;
    vector<double> strains;
};

// encapsulates all the data necessary to cache a simulation (sim data, synth pop, rngs)
class SimulationCache {
public:
    SimulationCache() {
        community     = nullptr;
        sim_ledger    = nullptr;
        rng           = nullptr;
        reporting_rng = nullptr;
    };

    SimulationCache(Community* o_community, SimulationLedger* o_sim_ledger, gsl_rng* o_rng, gsl_rng* o_reporting_rng) {
        community     = new Community(*o_community);
        sim_ledger    = new SimulationLedger(*o_sim_ledger);
        rng           = gsl_rng_clone(o_rng);
        reporting_rng = gsl_rng_clone(o_reporting_rng);
    };

    ~SimulationCache() {
        delete community;
        delete sim_ledger;
        gsl_rng_free(rng);
        gsl_rng_free(reporting_rng);
    };

    Community* community;
    SimulationLedger* sim_ledger;
    gsl_rng* rng;
    gsl_rng* reporting_rng;
};

// helper function to call simvis.R when needed
void gen_simvis(vector<string> &plot_log_buffer) {
    for (size_t i = 1; i < plot_log_buffer.size(); ++i) {
        //plot_log_buffer[i] = plot_log_buffer[i] + "," + to_string(Rt[i-1].second);
        if (i >= (plot_log_buffer.size() - 14)) {
            plot_log_buffer[i] = plot_log_buffer[i] + ",0";// + to_string(Rt_ma[i-1]);
        } else {
            plot_log_buffer[i] = plot_log_buffer[i];// + to_string(Rt_ma[i-1]);
        }
    }
    bool overwrite = true;
    write_daily_buffer(plot_log_buffer, "42", "plot_log.csv", overwrite);
    //int retval = system("Rscript fitvis.R");
    int retval = system("Rscript simvis.R");
    if (retval == -1) { cerr << "System call to `Rscript simvis.R` failed\n"; }
}


// helper function to handle score the fit in a tuning window
// NOTE: USING DEATHS BY DAY OF DEATH NOT DAY OF REPORTING
double score_fit(const Parameters* par, const Community* community, const size_t sim_day, vector<size_t> sim_rcases_szt, const map<size_t, vector<int>> &emp_data_map) {
    const int fit_window_size = (par->num_preview_windows + 1) * par->tuning_window;    // size of the tuning window
    const int window_start = ((int) sim_day + 1) - fit_window_size;                     // sim_day at the start of the window

    // create an interpolated curve between the error weights within the window
    vector<TimeSeriesAnchorPoint> fit_and_prev_anchors;
    if (window_start == 0) {
        fit_and_prev_anchors = {{Date::to_ymd(window_start,                                      par), 0.0},
                                {Date::to_ymd((window_start - 1) + par->tuning_window,           par), 1.0},
                                {Date::to_ymd((window_start - 1) + (2 * par->tuning_window),     par), 0.5},
//                                {Date::to_ymd((window_start - 1) + (2 * par->tuning_window) + 1, par), 0.5},
                                {Date::to_ymd(sim_day,                                           par), 0.0}};
    } else {
        fit_and_prev_anchors = {{Date::to_ymd(0,                                                 par), 0.0},
                                {Date::to_ymd(window_start,                                      par), 0.0},
                                {Date::to_ymd((window_start - 1) + par->tuning_window,           par), 1.0},
                                {Date::to_ymd((window_start - 1) + (2 * par->tuning_window),     par), 0.5},
//                                {Date::to_ymd((window_start - 1) + (2 * par->tuning_window) + 1, par), 0.5},
                                {Date::to_ymd(sim_day,                                           par), 0.0}};
    }
//    if (window_start == 0) { fit_and_prev_anchors.erase(fit_and_prev_anchors.begin()); }

    vector<double> fit_weights = Date::linInterpolateTimeSeries(fit_and_prev_anchors, par->startJulianYear, par->startDayOfYear);

    // cuts reported cases from the parsed empirical data
    vector<double> emp_rcases(sim_day + 1, 0.0); // default to 0 if data is missing, for purposes of calculating cumulative sum
    for (size_t day = 0; day <= sim_day; ++day) {
        if (emp_data_map.count(day)) {
            emp_rcases[day] = emp_data_map.at(day).at(0);
        } else {
            fit_weights[day] = 0.0;             // if there is no empirical data for this day, change the error weight to 0 (ignore the error on this day)
        }
    }

    const double sim_p10k = 1e4/community->getNumPeople();
    const double fl_pop   = 21538187;
    const double fl_p10k  = 1e4/fl_pop;

    for (auto& val: emp_rcases) { val *= fl_p10k; }     // normalize emp rcases to per 10000 people

    vector<double> sim_rcases(sim_rcases_szt.size());
    for (size_t i = 0; i < sim_rcases.size(); ++i) { sim_rcases[i] = sim_rcases_szt[i] * sim_p10k; }    // normalize sim rcases to per 10000 people

    // calculate cumulative sum of emp and sim rcases
    vector<double> emp_rcases_cumul(emp_rcases.size());
    vector<double> sim_rcases_cumul(sim_rcases.size());
    partial_sum(emp_rcases.begin(), emp_rcases.end(), emp_rcases_cumul.begin());
    partial_sum(sim_rcases.begin(), sim_rcases.end(), sim_rcases_cumul.begin());

    double distance = 0.0;
    double error = 0.0;
    double abs_distance = 0.0;
    double abs_error = 0.0;

    // for each day from the beginning of the sim to now, calculate the error (sim rcases - emp rcases) and the distance (error * weight)
    for (size_t d = 0; d <= sim_day; ++d) {
        //string date = Date::to_ymd(d, par);
        const double daily_crcase_error    = sim_rcases_cumul[d] - emp_rcases_cumul[d];
        const double daily_crcase_distance = daily_crcase_error * fit_weights[d];
        error    += daily_crcase_error;
        distance += daily_crcase_distance;

        abs_error    += abs(daily_crcase_error);
        abs_distance += abs(daily_crcase_distance);
    }

    // handles normalizing distance based on avg rep case incidence (offset by given number of days)
    const size_t fit_norm_offset = par->tuning_window * par->num_preview_windows;   // offset (in days) for when to calculate avg rep case incidence
    size_t num_days_to_avg = 3; // might consider 7 (one week)
    size_t epi_3day_window_start = (int) sim_day - fit_norm_offset - (num_days_to_avg - 1)/2;   // day to start averaging at
    vector<size_t> rcases_to_avg(emp_rcases.begin() + epi_3day_window_start, emp_rcases.begin() + epi_3day_window_start + num_days_to_avg);
    const double mean_rcases  = mean(rcases_to_avg);

    vector<double> single_rcase_p10k = {0, 0, (1 * fl_p10k)};   // for use if there are no rep cases to average (needed to prevent division by 0 later on)

    double avg_3day_rcases_w_offset = mean_rcases == 0 ? mean(single_rcase_p10k) : mean_rcases;

    const double normed_distance = distance/avg_3day_rcases_w_offset;
    return normed_distance;
}


// helps to keep track of the range the binary search can still act on
class Range {
public:
    Range() {};
    Range(double _min, double _max) : min(_min), max(_max) {};
//    Range(initializer_list<double> vals) : min(*(vals.begin())), max(*(vals.begin()+1)) {};

    void set_range(double _min, double _max) {
        min = _min;
        max = _max;
    }

    double min;
    double max;

    double best_value;
    double best_distance;
};


// encapsulates data and helper functions used to coordinate the auto tuning system
class BehaviorAutoTuner {
  public:
    BehaviorAutoTuner() {};
    ~BehaviorAutoTuner() { delete bin_search_range; };

    void init_defaults() {
        bin_search_range = new Range(0.0, 1.0);
        cur_anchor_val = 0.0;
        window_already_scored = false;

        tuning_window_ct = 1;

        manual_control = true;
        slow_auto = true;
        usr_choice = '\0';

        recache = false;

        fit_threshold = 20.0;

        output_buffer.str("");
    }

    void user_sim_setup() {
        cerr << "DO YOU WANT TO SIMULATE MANUALLY? (y/n) ";
        cin >> usr_choice;
        manual_control = (usr_choice == 'y') ? true :
                         (usr_choice == 'n') ? false : true;

        if (not manual_control) {
            cerr << "DO YOU WANT THE AUTO FITTING TO WAIT FOR KEYPRESSES BEFORE CONTINUING? (y/n) ";
            cin >> usr_choice;
            slow_auto = (usr_choice == 'y') ? true :
                        (usr_choice == 'n') ? false : true;
        }
    }

    double bin_search_range_min() { return bin_search_range->min; }
    double bin_search_range_max() { return bin_search_range->max; }

    void reset_bin_search_range() {
        bin_search_range->min = 0.0;
        bin_search_range->max = 1.0;
        window_already_scored = false;
    }

    void update_best_distance(double new_distance) {
        if (window_already_scored) {
            if (abs(new_distance) < abs(best_distance)) {
                best_distance = new_distance;
                best_anchor_val = cur_anchor_val;
            }
        } else {
            best_distance = new_distance;
            best_anchor_val = cur_anchor_val;
            window_already_scored = true;
        }
    }

    void clear_output_buffer() { output_buffer.str(""); }

    void print_header() {
        output_buffer << right << "  day          window  emp data%  sum dist   cur val  new best              search range   new val";
        cerr << right << output_buffer.str() << endl;
        clear_output_buffer();
    }

    Range* bin_search_range;        // points to the range used for this tuner's binary searches
    double cur_anchor_val;          // keeps track of the current anchor val being tested

    double best_anchor_val;
    double best_distance;
    bool window_already_scored;

    size_t tuning_window_ct;        // keeps a count of how many windows have been tuned as the simulation progresses

    bool manual_control;            // allow user input?
    bool slow_auto;                 // wait for 'enter' keypress after each automated decision
    char usr_choice;                // stores user inputs

    bool recache;                   // controls when the simulation cahce needs to be updated

    double fit_threshold;           // set value to determine if a window's fit is "good"

    map<size_t, vector<int>> emp_data;

    stringstream output_buffer;     // buffer to store data to output to the screen during auto tuning
};


// helper function to conduct modified binary search
double bin_search_anchor(BehaviorAutoTuner* tuner, double distance) {
    // grab the current search range and anchor value from the tuner
    Range* range = tuner->bin_search_range;
    const double cur_val = tuner->cur_anchor_val;
//    if (abs(range.best_distance) > abs(distance) ) {
//        range.best_distance = distance;
//        range.best_value    = cur_val;
//    }

    // prevents infinite searching; if a proposed anchor is less than MIN_ADJ different than cur_val, we stop the search
    const double MIN_ADJ = 0.01;
    double new_val;

    // we should never call for a search in these cases (this should have been caught in the calling scope)
    assert(cur_val <= range->max);
    assert(cur_val >= range->min);
    assert(range->max >= range->min);

    if (distance > 0) {
        range->set_range(cur_val, range->max);  // if distance is positive, the current val becomes the new search min
    } else {
        range->set_range(range->min, cur_val);  // if distance is negative, the current val becomes the new search max
    }

    const double adj = (range->max - range->min) / 2.0;     // calculates the adjustment to make on the current val
    // if adj is less than MIN_ADJ, we use the range max/min depending on the sign of the distance
    if (adj < MIN_ADJ) {
        if (distance > 0) {
            new_val = range->max;
            range->min = range->max;
        } else {
            new_val = range->min;
            range->max = range->min;
        }
    } else {
        // if the adj is greater than MIN_ADJ, we use the mid-point of the updated range as the new val
        new_val = range->min + adj;
    }
    return new_val;
}


// if auto tuning, this is the first step (create a new tuner, initialize defaults, ask user if manual control is desired)
BehaviorAutoTuner* initialize_behavior_autotuning(const Parameters* par, string emp_data_filename) {
    BehaviorAutoTuner* tuner = new BehaviorAutoTuner();
    tuner->init_defaults();
    tuner->user_sim_setup();

    tuner->emp_data = parse_emp_data_file(par, emp_data_filename);

    return tuner;
}


// necessary helper function since the first tuning window requires two anchor vals
void first_tuning_window_setup(const Parameters* par, Community* community, BehaviorAutoTuner* tuner, vector<TimeSeriesAnchorPoint> &social_distancing_anchors) {
    social_distancing_anchors.clear();
    community->_clearSocialDistancingTimedIntervention();

    double val1, val2;
    if (tuner->manual_control) {
        cerr << "ENTER SOC_CONTACT PARAMS FOR FIRST AND SECOND ANCHOR POINTS: ";
        cin >> val1 >> val2;
    } else {
        val1 = 0.0;//0.25;
        val2 = val1;
        tuner->cur_anchor_val = val1;

        tuner->print_header();
    }

    social_distancing_anchors.emplace_back(Date::to_ymd(0, par), val1);
    social_distancing_anchors.emplace_back(Date::to_ymd(par->tuning_window - 1, par), val2);
    community->setSocialDistancingTimedIntervention(social_distancing_anchors);
}


void overwrite_sim_cache(SimulationCache* &sim_cache, Community* community, SimulationLedger* ledger, BehaviorAutoTuner* tuner) {
    if (sim_cache) { delete sim_cache; }
    sim_cache = new SimulationCache(community, ledger, RNG, REPORTING_RNG);
    tuner->tuning_window_ct++;
    tuner->recache = false;
}


// helper function to handle main processing of a tuning window
void process_behavior_fit(int fit_is_good, double fit_distance, const Parameters* par, BehaviorAutoTuner* tuner, vector<TimeSeriesAnchorPoint> &social_distancing_anchors) {
    // add more data to the output
    string best_str = tuner->cur_anchor_val == tuner->best_anchor_val ? "  *       " : "          ";
    tuner->output_buffer << right
                         << setprecision(2) << scientific << setw(10) << fit_distance
                         << setprecision(6) << defaultfloat << setw(10) << tuner->cur_anchor_val
                         << best_str
                         << setw(3) << "[" << setw(10) << tuner->bin_search_range->min << ", " << setw(10) << tuner->bin_search_range->max << "]";

    if(fit_is_good) { // calling scope deemed window's fit as "good"
        tuner->recache = true;  // because fit is good, we can overwrite cache

        double val1;
        if (tuner->manual_control) {
            cerr << "FIT IS GOOD" << endl << "Enter next soc_contact step: ";
            cin >> val1;
        } else {
            val1 = tuner->cur_anchor_val;   // always start the next window with the val that was selected
            tuner->reset_bin_search_range();
            if (tuner->slow_auto) { cin.ignore(); }

            tuner->output_buffer << right << setw(10) << tuner->cur_anchor_val;
            cerr << right << tuner->output_buffer.str() << endl;
            tuner->clear_output_buffer();
            cerr << endl << endl;
            tuner->print_header();
        }

        social_distancing_anchors.emplace_back(Date::to_ymd((par->tuning_window * (tuner->tuning_window_ct + 1)) - 1, par), val1);
    } else { // calling scope deemed window's fit as "not good"
        // if using the automated system, conduct modified binary search for new anchor val
        if (not tuner->manual_control) { tuner->cur_anchor_val = bin_search_anchor(tuner, fit_distance); }

        if (tuner->tuning_window_ct == 1) { // tuning the first window
            double val1, val2;
            if (tuner->manual_control) {
                cerr << "FIT IS NOT GOOD" << endl << "Enter soc_contact params for first and second anchor points: ";
                cin >> val1 >> val2;
            } else {
                tuner->output_buffer << right << setw(10) << tuner->cur_anchor_val;
                val1 = tuner->cur_anchor_val;
                val2 = val1;
            }

            social_distancing_anchors[0] = {Date::to_ymd(0, par), val1};
            social_distancing_anchors[1] = {Date::to_ymd(par->tuning_window - 1, par), val2};
        } else { // past the first tuning window
            double val;
            if (tuner->manual_control) {
                cerr << "Enter new value for this window's soc_contact param: ";
                cin >> val;
            } else {
                tuner->output_buffer << right << setw(10) << tuner->cur_anchor_val;
                val = tuner->cur_anchor_val;
            }

            social_distancing_anchors.back() = {Date::to_ymd((par->tuning_window * tuner->tuning_window_ct) - 1, par), val};
        }
        //  this is what "waits" for the 'enter' keypress
        if (not tuner->manual_control and tuner->slow_auto) { cin.ignore(); }

        // output tuning data to the screen
        cerr << right << tuner->output_buffer.str() << endl;
        tuner->clear_output_buffer();
     }
}


void restore_from_cache(Community* &community, Date* &date, SimulationCache* sim_cache, SimulationLedger* &ledger, vector<TimeSeriesAnchorPoint> social_distancing_anchors) {
    if (ledger) { delete ledger; }
    ledger = new SimulationLedger(*(sim_cache->sim_ledger));

    if (community) { delete community; }
    community = new Community(*(sim_cache->community));
    date      = community->get_date();
    community->setSocialDistancingTimedIntervention(social_distancing_anchors);

    if (sim_cache->rng)           { gsl_rng_memcpy(RNG, sim_cache->rng); }
    if (sim_cache->reporting_rng) { gsl_rng_memcpy(REPORTING_RNG, sim_cache->reporting_rng); }
}


// main helper function to control the behavior auto tuning system
void behavior_autotuning(const Parameters* par, Community* &community, Date* &date, SimulationLedger* &ledger, BehaviorAutoTuner* tuner, SimulationCache* &sim_cache, vector<TimeSeriesAnchorPoint> &social_distancing_anchors) {
    const size_t day = date->day();

    if (tuner->recache and (day + 1) == tuner->tuning_window_ct * par->tuning_window) {
        // we are at the end of tuning window that was deemed "good"
        overwrite_sim_cache(sim_cache, community, ledger, tuner);
    } else if ((day + 1) == (tuner->tuning_window_ct + par->num_preview_windows) * par->tuning_window) {
        // calculates the proportion of days in the tuning window that has empirical data
        size_t window_start_sim_day = (day + 1) - ((par->num_preview_windows + 1) * par->tuning_window);
        size_t num_days_w_emp_data = 0;
        for (size_t d = window_start_sim_day; d <= day; ++d) { if (tuner->emp_data.count(d)) { ++num_days_w_emp_data; } }
        double prop_days_w_emp_data = (((double) num_days_w_emp_data) / (day - window_start_sim_day + 1));

        // begin adding tuning data to print to screen
        tuner->output_buffer << right
                             << setw(5) << day
                             << setw(3) << "[" << setw(5) << window_start_sim_day << ", " << setw(5) << day << "]"
                             << setw(10) << prop_days_w_emp_data * 100 << "%";

        // calculate normalized distance for this tuning window
        double fit_distance = score_fit(par, community, day, community->getNumDetectedCasesReport(), tuner->emp_data);
//cerr << "current anchor, score: " << tuner->cur_anchor_val << ", " << fit_distance << endl;
        // keep track of the anchor val with the smallest distance for this window
        tuner->update_best_distance(fit_distance);

        gen_simvis(ledger->plot_log_buffer);

        int fit_is_good = 0;
        if (tuner->manual_control) {
            cerr << "IS THE FIT GOOD? ";
            cin >> fit_is_good;
        } else {
            //fit_is_good = abs(fit_distance) < tuner->fit_threshold
            fit_is_good = abs(fit_distance) < (tuner->fit_threshold * prop_days_w_emp_data)                 // is the abs(distance) less than the threshold after adjusting for the presence of emp data
                          or (fit_distance > 0 and tuner->cur_anchor_val == tuner->bin_search_range_max())  // is the distance pos and the cur val is the max of the search range (can't search more)
                          or (fit_distance < 0 and tuner->cur_anchor_val == tuner->bin_search_range_min())  // is the distance neg and the cur val is the min of the search range (can't search more)
                          or (tuner->tuning_window_ct > 1 and prop_days_w_emp_data == 0)
                          or (tuner->bin_search_range_max() == tuner->bin_search_range_min());              // are the search range min and max equal (can't search more)

            if (fit_is_good) {
                // if the fit is deemed "good" ensure that we will use the best val found
                tuner->cur_anchor_val = tuner->best_anchor_val;
                fit_distance = tuner->best_distance;
cerr << "Decision made." << endl;
//cerr << "Using anchor val: " << tuner->cur_anchor_val << endl;
                if (tuner->tuning_window_ct > 1 and prop_days_w_emp_data == 0.25) {
                    tuner->cur_anchor_val = min_element(community->getTimedIntervention(SOCIAL_DISTANCING));
                }
            }
        }

        process_behavior_fit(fit_is_good, fit_distance, par, tuner, social_distancing_anchors);
        restore_from_cache(community, date, sim_cache, ledger, social_distancing_anchors);
    }
}


// if using previous tuned values, parse that file and save the behavior vals
void init_behavioral_vals_from_file(const Parameters* par, Community* community) {
    vector< vector<string> > tuning_data = read_2D_vector_file(par->autotuning_dataset, ',');
    vector<double> behavior_vals;
    vector<pair<int,double>> vals_to_interpolate;
    bool interpolating = false;

    for (vector<string> &v : tuning_data) {
        if (v[0] == "date") { continue; }
        if (v[0] == "interpolate") { interpolating = true; continue; }

        if (interpolating) {
            int sim_day = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, v[0]);
            double val  = stod(v[3]);
            vals_to_interpolate.push_back({sim_day, val});
        } else {
            double behavior = stod(v[3]);
            behavior_vals.push_back(behavior);
        }
    }

    if (interpolating) {
        if (vals_to_interpolate.size() > 1) {
            for (size_t i = 1; i < vals_to_interpolate.size(); ++i) {
                int time_delta   = vals_to_interpolate[i].first - vals_to_interpolate[i - 1].first;
                double start_val = vals_to_interpolate[i - 1].second;
                double end_val   = vals_to_interpolate[i].second;
                vector<double> interpolated_behavior = Date::linInterpolate(start_val, end_val, time_delta);
                behavior_vals.insert(behavior_vals.end(), interpolated_behavior.begin(), interpolated_behavior.end());
            }
        } else {
            cerr << "ERROR: to interpolate PPB, more than one anchor point is needed in " << par->autotuning_dataset << endl;
            exit(-2);
        }
    }

    community->setSocialDistancingTimedIntervention(behavior_vals);
}


vector<string> simulate_epidemic(const Parameters* par, Community* &community, const string process_id, const vector<string> mutant_intro_dates) {
    SimulationLedger* ledger    = new SimulationLedger();
    SimulationCache* sim_cache  = nullptr;

    BehaviorAutoTuner* tuner    = nullptr;

    //bool window_fit_is_good = true;
    //vector<int> epi_sizes;
    //vector<string> daily_output_buffer;

    //map<string, vector<int> > periodic_incidence = construct_tally();
    ledger->periodic_incidence = construct_tally();
    //vector<int> periodic_prevalence(NUM_OF_PREVALENCE_REPORTING_TYPES, 0);
    ledger->periodic_prevalence = vector<int>(NUM_OF_PREVALENCE_REPORTING_TYPES, 0);
    vector<double> trailing_averages(par->runLength);
    const double pop_at_risk = min(community->getNumPeople(), par->numSurveilledPeople);

    // default behavioral anchor points
    vector<TimeSeriesAnchorPoint> social_distancing_anchors = {
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
        {"2021-10-01", 0.0}
    };
    community->setSocialDistancingTimedIntervention(social_distancing_anchors);

    //vector<string> plot_log_buffer = {"date,sd,seasonality,vocprev1,vocprev2,cinf,closed,rcase,rdeath,inf,rhosp,Rt"};
    ledger->plot_log_buffer = {"serial,date,sd,seasonality,vocprev1,vocprev2,vocprev3,cinf,closed,rcase,rdeath,inf,rhosp,VES,brkthruRatio,vaxInfs,unvaxInfs,hospInc,hospPrev,icuInc,icuPrev,vaxHosp,unvaxHosp,Rt"};
    //ledger->plot_log_buffer = {"date,sd,seasonality,vocprev1,vocprev2,cinf,closed,rcase,rdeath,inf,rhosp,Rt"};

    Date* date = community->get_date();
    ledger->strains = {50.0, 0.0, 0.0, 0.0}; // initially all WILDTYPE
    assert(ledger->strains.size() == NUM_OF_STRAIN_TYPES);

    if (par->behavioral_autotuning) {
        // create tuner and initialize first simulation cache
        tuner = initialize_behavior_autotuning(par, "rcasedeath-florida.csv");
        sim_cache = new SimulationCache(community, ledger, RNG, REPORTING_RNG);
        first_tuning_window_setup(par, community, tuner, social_distancing_anchors);
    } else if (not par->autotuning_dataset.empty()) {
        // filename provided for a dataset with behavioral values to use
        init_behavioral_vals_from_file(par, community);
    }

    for (; date->day() < (signed) par->runLength; date->increment()) {
if (*date == "2021-12-01") { gsl_rng_set(RNG, par->randomseed); }
        community->tick();

        if (par->behavioral_autotuning) { behavior_autotuning(par, community, date, ledger, tuner, sim_cache, social_distancing_anchors); }

        //update_vaccinations(par, community, date);
        //community->tick();
        const size_t sim_day = date->day();

        if ( mutant_intro_dates.size() ) {
            if (*date >= mutant_intro_dates[0] and *date < mutant_intro_dates[1]) {
                //const int time_since_intro = date->day() - Date::to_sim_day(par->julian_start_day, par->julian_start_year, mutant_intro_dates[0]);
                if (ledger->strains[WILDTYPE] > 1) {
                    ledger->strains[WILDTYPE]--;
                    ledger->strains[ALPHA]++;
                }
            } else if (*date >= mutant_intro_dates[1] and *date < mutant_intro_dates[2]) {
                if (ledger->strains[WILDTYPE] > 1) {
                    ledger->strains[WILDTYPE]--;
                    ledger->strains[DELTA]++;
                }
                if (ledger->strains[ALPHA] > 1) {
                    ledger->strains[ALPHA]--;
                    ledger->strains[DELTA]++;
                }
            } else if (*date >= mutant_intro_dates[2]) {
                for (int i = 0; i < 5; ++i) { // faster take-over of omicron
                    if (ledger->strains[WILDTYPE] > 1) {
                        ledger->strains[WILDTYPE]--;
                        ledger->strains[OMICRON]++;
                    }
                    if (ledger->strains[ALPHA] > 1) {
                        ledger->strains[ALPHA]--;
                        ledger->strains[OMICRON]++;
                    }
                    if (ledger->strains[DELTA] > 1) {
                        ledger->strains[DELTA]--;
                        ledger->strains[OMICRON]++;
                    }
                }
            }
        }

        seed_epidemic(par, community, date, ledger->strains);
        const vector<size_t> infections         = community->getNumNewlyInfected();
        const vector<size_t> all_reported_cases = community->getNumDetectedCasesReport();
        const size_t reported_cases             = all_reported_cases[sim_day];
        //const double trailing_avg               = calc_trailing_avg(all_reported_cases, sim_day, 7); // <= 7-day trailing average
        const vector<size_t> rhosp              = community->getNumDetectedHospitalizations();
        const vector<size_t> severe_prev        = community->getNumSeverePrev();
        const double cinf                       = accumulate(infections.begin(), infections.begin()+sim_day+1, 0.0);
        const double cAR                        = cinf/pop_at_risk; // cumulative attack rate (I hate this term)
        const vector<size_t> rdeaths            = community->getNumDetectedDeathsOnset();

        const vector<size_t> severe             = community->getNumNewlySevere();
        const double trailing_avg = trailing_averages[sim_day];

        const size_t rc_ct = accumulate(all_reported_cases.begin(), all_reported_cases.begin()+sim_day+1, 0);
        map<string, double> VE_data = community->calculate_daily_direct_VE();
        //vector<string> inf_by_loc_keys = {"home", "social", "work_staff", "patron", "school_staff", "student", "hcw", "patient", "ltcf_staff", "ltcf_resident"};
        //for (string key : inf_by_loc_keys) {
        //    cerr << "infLoc " << date->to_ymd() << ' ' << key << ' ' << community->getNumNewInfectionsByLoc(key)[sim_day] << endl;
        //}
        if (not par->behavioral_autotuning) {
            if (date->dayOfMonth()==1) cerr << "        rep sday        date  infinc  cAR     rcases  rcta7  crcases  rdeath  crdeath  sevprev   crhosp  closed  socdist  coverage\n";
            cerr << right
                << setw(11) << "NA" //process_id
                << setw(5)  << sim_day
                << setw(12) << date->to_ymd()
                << setw(8)  << infections[sim_day]
                << "  "     << setw(7) << setprecision(2) << left << cAR << right
                << setw(7)  << reported_cases
                << "  "     << setw(7) << setprecision(4) << left << trailing_avg << right
                << setw(7)  << rc_ct
                << setw(8)  << rdeaths[sim_day]
                << setw(9)  << accumulate(rdeaths.begin(), rdeaths.begin()+sim_day+1, 0)
                << setw(9)  << severe_prev[sim_day]
                << setw(9)  << accumulate(rhosp.begin(), rhosp.begin()+sim_day+1, 0)
                << setw(8)  << community->getTimedIntervention(NONESSENTIAL_BUSINESS_CLOSURE, sim_day)
                << "  "     << setprecision(2) << community->getTimedIntervention(SOCIAL_DISTANCING, sim_day)
                << "  "     << setprecision(2) << VE_data["coverage"]
                << "  "     << setprecision(2) << VE_data["dose_1"]
                << "  "     << setprecision(2) << VE_data["dose_2"]
                << "  "     << setprecision(2) << VE_data["dose_3"]
                //             << "  "     << setprecision(2) << (double) severe[sim_day] / reported_cases
                << endl;
        }

        //date,sd,seasonality,vocprev,cinf,closed,rcase,rdeath,Rt
        stringstream ss;
        ss << process_id << ","
           << date->to_string({"yyyy", "mm", "dd"}, "-") << ","
           << community->social_distancing(sim_day)/*par->timedInterventions.at(SOCIAL_DISTANCING).at(sim_day)*/ << ","
           << par->seasonality.at(date->julianDay()-1) << ","
           << (float) community->getNumNewInfections(ALPHA)[sim_day]/infections[sim_day] << ","
           << (float) community->getNumNewInfections(DELTA)[sim_day]/infections[sim_day] << ","
           << (float) community->getNumNewInfections(OMICRON)[sim_day]/infections[sim_day] << ","
           << cAR << ","
           << community->getTimedIntervention(NONESSENTIAL_BUSINESS_CLOSURE, sim_day)<< ","
           << reported_cases*1e4/pop_at_risk << ","
           << rdeaths[sim_day]*1e4/pop_at_risk << ","
           << infections[sim_day]*1e4/pop_at_risk << ","
           << rhosp[sim_day]*1e4/pop_at_risk << ","
           << VE_data["VES"] << ","
           << VE_data["breakthruRatio"] << ","
           << VE_data["vaxInfs"]*1e4/pop_at_risk << ","
           << VE_data["unvaxInfs"]*1e4/pop_at_risk << ","
           << community->getNumHospInc()[sim_day]*1e4/pop_at_risk << ","
           << community->getNumHospPrev()[sim_day]*1e4/pop_at_risk << ","
           << community->getNumIcuInc()[sim_day]*1e4/pop_at_risk << ","
           << community->getNumIcuPrev()[sim_day]*1e4/pop_at_risk << ","
           << VE_data["vaxHosp"]*1e4/pop_at_risk << ","
           << VE_data["unvaxHosp"]*1e4/pop_at_risk;
        ledger->plot_log_buffer.push_back(ss.str());

    }

    const vector<size_t> sim_reported_cases = community->getNumDetectedCasesReport();
    if (par->behavioral_autotuning) {
        ofstream ofs("autotuning_dataset.csv");
        ofs << "date,sim_rcase,emp_rcase,behavior" << endl;
        for (size_t i = 0; i < par->runLength; ++i) {
            int daily_emp_data = tuner->emp_data.count(i) ? tuner->emp_data.at(i)[0] : 0;
            ofs << Date::to_ymd(i, par) << "," << sim_reported_cases[i] << "," << daily_emp_data << "," << community->getTimedIntervention(SOCIAL_DISTANCING, i) << endl;
        }
        ofs.close();
    }

    if (sim_cache)  { delete sim_cache; }

    double cdeath_icu   = 0.0;
    double cdeath2      = 0.0;
    for (Person* p: community->getPeople()) {
        if (p->getNumNaturalInfections() and p->getInfection()->fatal()) {
            // this counts *all* deaths and icu admissions, including those that were scheduled
            // but didn't happen because the simulation ended during the infection
            cdeath2++;
            cdeath_icu += p->getInfection()->icu();
        }
    }

    const vector<size_t> infections         = community->getNumNewlyInfected();
    const vector<size_t> symptomatic        = community->getNumNewlySymptomatic();
    const vector<size_t> severe             = community->getNumNewlySevere();
    const vector<size_t> critical           = community->getNumNewlyCritical();
    const vector<size_t> dead               = community->getNumNewlyDead();
        const vector<size_t> all_reported_cases = community->getNumDetectedCasesReport();
//
    const double cinf  = accumulate(infections.begin(), infections.end(), 0.0);
    const double csymp = accumulate(symptomatic.begin(), symptomatic.end(), 0.0);
    const double csev  = accumulate(symptomatic.begin(), symptomatic.end(), 0.0);
//    const double ccrit = accumulate(symptomatic.begin(), symptomatic.end(), 0.0);
//    const double cdead = accumulate(symptomatic.begin(), symptomatic.end(), 0.0);
        const size_t rc_ct = accumulate(all_reported_cases.begin(), all_reported_cases.end(), 0);

    cerr << "symptomatic infections, total infections, asymptomatic fraction: " << csymp << ", " << cinf << ", " << 1.0 - (csymp/cinf) << endl;
    cerr << "icu deaths, total deaths, ratio: " << cdeath_icu << ", " << cdeath2 << ", " << cdeath_icu/cdeath2 << endl;
    cerr << "severe infections / all reported cases: " << (double) csev / rc_ct << endl;

cerr_vector(community->getTimedIntervention(SOCIAL_DISTANCING));

vector<size_t> offspringDistribution = community->generateOffspringDistribution();
cerr << endl << "OFFSPRING DISTR" << endl;
cerr_vector(offspringDistribution);
cerr << endl;

    if (RNG) { gsl_rng_free(RNG); }
    if (REPORTING_RNG) { gsl_rng_free(REPORTING_RNG); }

    vector<string> plot_log_buffer = ledger->plot_log_buffer;
    if (ledger) { delete ledger; }
    if (tuner)  { delete tuner; }


    return plot_log_buffer;
}


bool fileExists(const std::string& filename) {
    struct stat buf;
    return stat(filename.c_str(), &buf) != -1;
}


void write_daily_buffer(vector<string>& buffer, const string process_id, string filename, bool overwrite) {
    if (filename == "") {
        stringstream ss_filename;
        ss_filename << "daily_output." << process_id;
        filename = ss_filename.str();
    }

    string all_output;
    for (const auto &line : buffer) all_output += (line + "\n");

    if (fileExists(filename) and not overwrite) {
        cerr << "WARNING: Daily output file already exists: " << filename << endl << "WARNING: Aborting write.\n";
        return;
    }

    ofstream file;
    file.open(filename);

    if (file.is_open()) {  // TODO - add this check everywhere a file is opened
        file << all_output;
        file.close();
    } else {
        cerr << "ERROR: Could not open daily buffer file for output: " << filename << endl;
        exit(-842);
    }
}


void write_immunity_by_age_file(const Community* community, const int year, string filename) {
    if (filename == "") {
        stringstream ss_filename;
        ss_filename << "imm_vs_age.year" << year;
        filename = ss_filename.str();
    }
    // total count, immune
    vector< vector<int> > tally(NUM_AGE_CLASSES, vector<int>(2, 0));
    for (Person* p: community->getPeople()) {
        const int age = p->getAge();
        tally[age][0]++;
        const int numInfections = p->getNumNaturalInfections();
        if (numInfections > 0) tally[age][1]++;
    }

    ofstream file;
    file.open(filename);
    for (int a = 0; a<NUM_AGE_CLASSES; ++a) {
        file << a;
        for (int i: tally[a]) file << " " << i;
        file << endl;
    }
    file.close();
}

void write_immunity_file(const Community* community, const string label, string filename, int runLength) {
    if (filename == "") {
        stringstream ss_filename;
        ss_filename << "immunity." << label;
        filename = ss_filename.str();
    }
    ofstream file;
    file.open(filename);
    file << "pid age infection_history\n";
    for (Person* p: community->getPeople()) {
        file << p->getID() << " " << p->getAge();
        for (int k = 0; k<p->getNumNaturalInfections(); ++k) {
            const int infection_history = p->getInfectedTime(k) - runLength; // 0 is no infection; -1 means yesterday, -2 means 2 days ago ...
            file << " " << infection_history;
        }
        file << endl;
    }
    file.close();
}


void daily_detailed_output(Community* community, int t) {
    // print out infected people
    for (Person *p: community->getPeople()) {
        if (p->isInfected(t)) {
            // home location
            cout << t
                 << ",p,"
                 << p->getID() << ","
                 << p->getHomeLoc()->getID() << ","
                 << (p->isSymptomatic(t)?1:0) << ","
                 //<< (p->isWithdrawn(t)?1:0) << ","
                 << (p->isNewlyInfected(t)?1:0) << endl;
        }
    }
}


void generate_infection_db(const Community* community, const unsigned long int serial) {
    stringstream infection_filename;
    infection_filename << "./infection_history_" << serial << ".csv";

    stringstream offspring_filename;
    offspring_filename << "./offspring_infections_" << serial << ".csv";

    stringstream detection_filename;
    detection_filename << "./infection_detection_" << serial << ".csv";

    ofstream infection_file(infection_filename.str(), std::ios::trunc); // for the infection table
    ofstream offspring_file(offspring_filename.str(), std::ios::trunc); // for the offspring table
    ofstream detection_file(detection_filename.str(), std::ios::trunc); // for the detection table

    if (infection_file and offspring_file and detection_file) { cerr << "ALL FILES OPEN" << endl; }
    else { cerr << "FILES FAILED TO OPEN" << endl; }

    for (Person* p : community->getPeople()) {
        for (Infection* inf : p->getInfectionHistory()) {
            if (not inf) { continue; }
            int inf_place_id = inf->getInfectedPlace() ? inf->getInfectedPlace()->getID() : -1;
            int inf_by_id    = inf->getInfectedBy() ? inf->getInfectedBy()->getID() : -1;
            int inf_owner_id = inf->getInfectionOwner() ? inf->getInfectionOwner()->getID() : -1;

            infection_file << inf << ','
                           << inf_place_id << ','
                           << inf_by_id << ','
                           << inf_owner_id << ','
                           << inf->getInfectedTime() << ','
                           << inf->getInfectiousTime() << ','
                           << inf->getInfectiousEndTime() << ','
                           << inf->getSymptomTime() << ','
                           << inf->getSymptomEndTime() << ','
                           << inf->getSevereTime() << ','
                           << inf->getSevereEndTime() << ','
                           << inf->getHospitalizedTime() << ','
                           << inf->getCriticalTime() << ','
                           << inf->getCriticalEndTime() << ','
                           << inf->getIcuTime() << ','
                           << inf->getDeathTime() << ','
                           << inf->getStrain() << ','
                           << inf->getRelInfectiousness() << ','
                           << inf->getDetection() << ','
                           << inf->secondary_infection_tally() << endl;

            for (Infection* sec_inf : inf->get_infections_caused()) {
                offspring_file << inf << ',' << sec_inf << endl;
            }

            if (inf->getDetection()) {
                Detection* inf_det = inf->getDetection();
                detection_file << inf_det << ','
                               << inf << ','
                               << inf_det->detected_state << ','
                               << inf_det->reported_time << endl;
            }
        }
    }

    infection_file.close();
    offspring_file.close();
    detection_file.close();

    stringstream ss;
    ss << "sqlite3 infection_data_" << serial << ".sqlite '.read gen_infection_db.sql'";
    string cmd_str = ss.str();
    int retval = system(cmd_str.c_str());
    if (retval == -1) { cerr << "System call to `sqlite3 infection_data_0.db '.read gen_infection_db.sql'` failed\n"; }

    ss.str(string());
    ss << "rm " << infection_filename.str() << ' ' << offspring_filename.str() << ' ' << detection_filename.str();
    cmd_str = ss.str();
    retval = system(cmd_str.c_str());
    if (retval == -1) { cerr << "System call to delete infection csv files failed\n"; }
}
