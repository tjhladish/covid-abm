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
    if (!community->loadPopulation(par->populationFilename, par->comorbidityFilename)) {
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

// NOTE: USING DEATHS BY DAY OF DEATH NOT DAY OF REPORTING
vector<double> fitting_error(const Parameters* par, const Community* community, const size_t sim_day, vector<size_t> sim_rcases_szt, const map<size_t, vector<int>> &emp_data_map) {
    const int fit_window_size = (par->num_preview_windows + 1) * par->fitting_window;
    const int window_start = ((int) sim_day + 1) - fit_window_size;
    vector<TimeSeriesAnchorPoint> fit_and_prev_anchors;
    if (window_start == 0) {
        fit_and_prev_anchors = {{Date::to_ymd(window_start,                                       par), 0.0},
                                {Date::to_ymd((window_start - 1) + par->fitting_window,           par), 1.0},
                                {Date::to_ymd((window_start - 1) + (2 * par->fitting_window),     par), 0.5},
                                {Date::to_ymd((window_start - 1) + (2 * par->fitting_window) + 1, par), 0.5},
                                {Date::to_ymd(sim_day,                                            par), 0.5}};
    } else {
        fit_and_prev_anchors = {{Date::to_ymd(0,                                                  par), 0.0},
                                {Date::to_ymd(window_start,                                       par), 0.0},
                                {Date::to_ymd((window_start - 1) + par->fitting_window,           par), 1.0},
                                {Date::to_ymd((window_start - 1) + (2 * par->fitting_window),     par), 0.5},
                                {Date::to_ymd((window_start - 1) + (2 * par->fitting_window) + 1, par), 0.5},
                                {Date::to_ymd(sim_day,                                            par), 0.5}};
    }
//    if (window_start == 0) { fit_and_prev_anchors.erase(fit_and_prev_anchors.begin()); }

    vector<double> fit_weights = Date::linInterpolateTimeSeries(fit_and_prev_anchors, par->startJulianYear, par->startDayOfYear);

    vector<double> emp_rcases(sim_day + 1, 0.0); // default to 0 if data is missing, for purposes of calculating cumulative sum
    for (size_t day = 0; day <= sim_day; ++day) {
        if (emp_data_map.count(day)) {
            emp_rcases[day] = emp_data_map.at(day).at(0);
        } else {
            fit_weights[day] = 0.0;
        }
    }

    const double sim_p10k = 1e4/community->getNumPeople();
    const double fl_pop   = 21538187;
    const double fl_p10k  = 1e4/fl_pop;

    for (auto& val: emp_rcases) { val *= fl_p10k; }

    vector<double> sim_rcases(sim_rcases_szt.size());
    for (size_t i = 0; i < sim_rcases.size(); ++i) { sim_rcases[i] = sim_rcases_szt[i] * sim_p10k; }
    //for (auto& val: sim_rcases) { val *= sim_p10k; }

    vector<double> emp_rcases_cumul(emp_rcases.size());
    vector<double> sim_rcases_cumul(sim_rcases.size());
    partial_sum(emp_rcases.begin(), emp_rcases.end(), emp_rcases_cumul.begin());
    partial_sum(sim_rcases.begin(), sim_rcases.end(), sim_rcases_cumul.begin());

cerr << fixed << setprecision(3);
cerr << "DATE\t\t\tCR ERROR\tADJ ERROR (ADJ)" << endl;

    double distance = 0.0;
    double error = 0.0;
    double abs_distance = 0.0;
    double abs_error = 0.0;

    for (size_t d = 0; d <= sim_day; ++d) {
        //string date = Date::to_ymd(d, par);
        const double daily_crcase_error    = sim_rcases_cumul[d] - emp_rcases_cumul[d];
        const double daily_crcase_distance = daily_crcase_error * fit_weights[d];
        error    += daily_crcase_error;
        distance += daily_crcase_distance;

        abs_error    += abs(daily_crcase_error);
        abs_distance += abs(daily_crcase_distance);
    }

cerr << "\nTOTAL\t\t\t" << abs_error << " (" << error << ")" << endl;//<< "\t------\t" << abs_error << " (" << (abs_error/raw_abs_error) * 100 << "%)" << endl;
cerr << "ADJ TOTAL\t\t" << abs_distance << " (" << distance << ")" << endl;//<< "\t------\t" << abs_error << " (" << (abs_error/raw_abs_error) * 100 << "%)" << endl;

    const size_t fit_norm_offset = par->fitting_window * par->num_preview_windows;
    size_t num_days_to_avg = 3; // might consider 7 (one week)
    size_t epi_3day_window_start = (int) sim_day - fit_norm_offset - (num_days_to_avg - 1)/2;
    vector<size_t> rcases_to_avg(emp_rcases.begin() + epi_3day_window_start, emp_rcases.begin() + epi_3day_window_start + num_days_to_avg);
    const double mean_rcases  = mean(rcases_to_avg);

    vector<double> single_rcase_p10k = {0, 0, (1 * fl_p10k)};

    double avg_3day_rcases_w_offset = mean_rcases == 0 ? mean(single_rcase_p10k) : mean_rcases;
cerr << "AVG CENTERED ON DAY " << Date::to_ymd(sim_day - fit_norm_offset, par) << endl;
cerr << "AVG RCASES: " << avg_3day_rcases_w_offset << endl;
//cerr << "RCASES AVERAGED: ";
//for (size_t val : tmp_3day_rcases_p10k) { cerr << val << " "; }
//cerr << endl;
if (avg_3day_rcases_w_offset > 0) { cerr << "NORMALIZED ERROR: " << abs_distance/avg_3day_rcases_w_offset << " (" << distance/avg_3day_rcases_w_offset << ")" << endl; }
cerr << defaultfloat;

    return {abs_distance/avg_3day_rcases_w_offset, distance/avg_3day_rcases_w_offset};
}

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
};

double bin_search_anchor(Range* range, double cur_val, double distance) {
    const double MIN_ADJ = 0.01;
    double new_val;
cerr << "CUR BIN SEARCH RANGE: [" << range->min << ", " << range->max << "]" << endl;
cerr << "CUR VAL, fit distance:              " << cur_val << ", " << distance << endl;
    assert(cur_val <= range->max);
    assert(cur_val >= range->min);
    assert(range->max >= range->min);

    if (distance > 0) {
        range->set_range(cur_val, range->max);
    } else {
        range->set_range(range->min, cur_val);
    }

    const double adj = (range->max - range->min) / 2.0;
    if (adj < MIN_ADJ) {
        if (distance > 0) {
            new_val = range->max;
            range->min = range->max;
        } else {
            new_val = range->min;
            range->max = range->min;
        }
    } else {
        new_val = range->min + adj;
    }
cerr << "ADJ, NEW VAL:         " << adj << ", " << new_val << endl;
cerr << "NEW BIN SEARCH RANGE: [" << range->min << ", " << range->max << "]" << endl;
    return new_val;
}


class BehaviorAutoTuner {
  public:
    BehaviorAutoTuner() {};
    ~BehaviorAutoTuner() { delete bin_search_range; };

    void init_defaults() {
        bin_search_range = new Range(0.0, 1.0);
        cur_anchor_val = 0.0;

        fitting_window_ct = 1;

        manual_control = true;
        slow_auto = true;
        usr_choice = '\0';

        recache = false;

        fit_threshold = 20.0;
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
    }

    Range* bin_search_range;
    double cur_anchor_val;

    size_t fitting_window_ct;

    bool manual_control;
    bool slow_auto;
    char usr_choice;

    bool recache;

    double fit_threshold;

    map<size_t, vector<int>> emp_data;
};

BehaviorAutoTuner* initialize_behavior_auto_tuning(const Parameters* par, string emp_data_filename) {
    BehaviorAutoTuner* tuner = new BehaviorAutoTuner();
    tuner->init_defaults();
    tuner->user_sim_setup();

    tuner->emp_data = parse_emp_data_file(par, emp_data_filename);

    return tuner;
}

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
        tuner->cur_anchor_val  = val1;
    }
    social_distancing_anchors.emplace_back(Date::to_ymd(0, par), val1);
    social_distancing_anchors.emplace_back(Date::to_ymd(par->fitting_window - 1, par), val2);
    community->setSocialDistancingTimedIntervention(social_distancing_anchors);
}

void overwrite_sim_cache(SimulationCache* &sim_cache, Community* community, SimulationLedger* ledger, BehaviorAutoTuner* tuner) {
    if (sim_cache) { delete sim_cache; }
    sim_cache = new SimulationCache(community, ledger, RNG, REPORTING_RNG);
    tuner->fitting_window_ct++;
    tuner->recache = false;
}

void process_behavior_fit(int fit_is_good, vector<double> fit_error, const Parameters* par, BehaviorAutoTuner* tuner, vector<TimeSeriesAnchorPoint> &social_distancing_anchors) {
    if(fit_is_good) {
cerr << "FIT IS GOOD" << endl;
        tuner->recache = true;

        double val1;
        if (tuner->manual_control) {
            cerr << "Enter next soc_contact step: ";
            cin >> val1;
        } else {
            cerr << tuner->cur_anchor_val << endl;
            val1 = tuner->cur_anchor_val;
            tuner->reset_bin_search_range();
            if (tuner->slow_auto) { cin.ignore(); }
        }

        social_distancing_anchors.emplace_back(Date::to_ymd((par->fitting_window * (tuner->fitting_window_ct + 1)) - 1, par), val1);
    } else {
cerr << "FIT IS NOT GOOD" << endl;
        if (not tuner->manual_control) { tuner->cur_anchor_val = bin_search_anchor(tuner->bin_search_range, tuner->cur_anchor_val, fit_error[1]); }

        if (tuner->fitting_window_ct == 1) {
            double val1, val2;
            if (tuner->manual_control) {
                cerr << "Enter soc_contact params for first and second anchor points: ";
                cin >> val1 >> val2;
            } else {
                cerr << tuner->cur_anchor_val << endl;
                val1 = tuner->cur_anchor_val;
                val2 = val1;
            }

            social_distancing_anchors[0] = {Date::to_ymd(0, par), val1};
            social_distancing_anchors[1] = {Date::to_ymd(par->fitting_window - 1, par), val2};
        } else {
            double val;
            if (tuner->manual_control) {
                cerr << "Enter new value for this window's soc_contact param: ";
                cin >> val;
            } else {
                cerr << tuner->cur_anchor_val << endl;
                val = tuner->cur_anchor_val;
            }

            social_distancing_anchors.back() = {Date::to_ymd((par->fitting_window * tuner->fitting_window_ct) - 1, par), val};
        }
        if (not tuner->manual_control and tuner->slow_auto) { cin.ignore(); }
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

void behavior_auto_tuning(const Parameters* par, Community* &community, Date* &date, SimulationLedger* &ledger, BehaviorAutoTuner* tuner, SimulationCache* &sim_cache, vector<TimeSeriesAnchorPoint> &social_distancing_anchors) {
    const size_t day = date->day();
    if (tuner->recache and (day + 1) == tuner->fitting_window_ct * par->fitting_window) {
        overwrite_sim_cache(sim_cache, community, ledger, tuner);
cerr << "RECACHING ON DAY " << day + 1 << endl;
    } else if ((day + 1) == (tuner->fitting_window_ct + par->num_preview_windows) * par->fitting_window) {
size_t window_start_sim_day = (day + 1) - ((par->num_preview_windows + 1) * par->fitting_window);
cerr << "FITTING WINDOW: " << window_start_sim_day << " TO " << window_start_sim_day + par->fitting_window - 1 << endl;
cerr << "PREVIEW WINDOW: " << window_start_sim_day + par->fitting_window  << " TO " << day << endl;
        vector<double> fit_error = fitting_error(par, community, day, community->getNumDetectedCasesReport(), tuner->emp_data);

        gen_simvis(ledger->plot_log_buffer);

        int tmp_cnt = 0;
        for (auto anchor: social_distancing_anchors) {
            if (tmp_cnt == 0) { cerr << anchor.value << " "; tmp_cnt++; }
            else { cerr << anchor.value << " 1 "; }
        }
        cerr << endl;
        cerr << "IS THE FIT GOOD? ";

        int fit_is_good = 0;
        if (tuner->manual_control) {
            cin >> fit_is_good;
        } else {
            cerr << endl;
            fit_is_good = abs(fit_error[1]) < tuner->fit_threshold
                          or (fit_error[1] > 0 and tuner->cur_anchor_val == tuner->bin_search_range_max())
                          or (fit_error[1] < 0 and tuner->cur_anchor_val == tuner->bin_search_range_min())
                          or (tuner->bin_search_range_max() == tuner->bin_search_range_min());
        }

        process_behavior_fit(fit_is_good, fit_error, par, tuner, social_distancing_anchors);
        restore_from_cache(community, date, sim_cache, ledger, social_distancing_anchors);
    }
}

vector<string> simulate_epidemic(const Parameters* par, Community* &community, const string /*process_id*/, const vector<string> mutant_intro_dates) {//,
                                 //map<size_t, TimeSeriesAnchorPoint> &social_contact_map) {
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
    ledger->plot_log_buffer = {"date,sd,seasonality,vocprev1,vocprev2,cinf,closed,rcase,rdeath,inf,rhosp,Rt"};

    Date* date = community->get_date();
    ledger->strains = {50.0, 0.0, 0.0}; // initially all WILDTYPE
    assert(ledger->strains.size() == NUM_OF_STRAIN_TYPES);

    if (par->auto_fitting) { // change to par->behavioral_auto_tuning
        tuner = initialize_behavior_auto_tuning(par, "rcasedeath-florida.csv");
        // inital sim cache
        sim_cache = new SimulationCache(community, ledger, RNG, REPORTING_RNG);

        first_tuning_window_setup(par, community, tuner, social_distancing_anchors);
    }

    for (; date->day() < (signed) par->runLength; date->increment()) {
        community->tick();

        if (par->auto_fitting) { behavior_auto_tuning(par, community, date, ledger, tuner, sim_cache, social_distancing_anchors); }

        //update_vaccinations(par, community, date);
        //community->tick();
        const size_t sim_day = date->day();

        if ( mutant_intro_dates.size() ) {
            if (*date >= mutant_intro_dates[0] and *date < mutant_intro_dates[1]) {
                //const int time_since_intro = date->day() - Date::to_sim_day(par->julian_start_day, par->julian_start_year, mutant_intro_dates[0]);
                if (ledger->strains[WILDTYPE] > 1) {
                    ledger->strains[WILDTYPE]--;
                    ledger->strains[B_1_1_7]++;
                }
            } else if (*date >= mutant_intro_dates[1]) {
                if (ledger->strains[WILDTYPE] > 1) {
                    ledger->strains[WILDTYPE]--;
                    ledger->strains[B_1_617_2]++;
                }
                if (ledger->strains[B_1_1_7] > 1) {
                    ledger->strains[B_1_1_7]--;
                    ledger->strains[B_1_617_2]++;
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
        //const vector<size_t> rdeaths            = community->getNumDetectedDeaths();
        const vector<size_t> rdeaths            = community->getNumDetectedDeathsOnset();


//        const size_t rc_ct = accumulate(all_reported_cases.begin(), all_reported_cases.begin()+sim_day+1, 0);
//        if (date->dayOfMonth()==1) cerr << "        rep sday        date  infinc  cAR     rcases  rcta7  crcases  rdeath  crdeath  sevprev   crhosp  closed  socdist\n";
//        cerr << right
//             << setw(11) << process_id
//             << setw(5)  << sim_day
//             << setw(12) << date->to_ymd()
//             << setw(8)  << infections[sim_day]
//             << "  "     << setw(7) << setprecision(2) << left << cAR << right
//             << setw(7)  << reported_cases
//             << "  "     << setw(7) << setprecision(4) << left << trailing_avg << right
//             << setw(7)  << rc_ct
//             << setw(8)  << rdeaths[sim_day]
//             << setw(9)  << accumulate(rdeaths.begin(), rdeaths.begin()+sim_day+1, 0)
//             << setw(9)  << severe_prev[sim_day]
//             << setw(9)  << accumulate(rhosp.begin(), rhosp.begin()+sim_day+1, 0)
//             << setw(8)  << community->getTimedIntervention(NONESSENTIAL_BUSINESS_CLOSURE, sim_day)
//             << "  "     << setprecision(2) << community->social_distancing(sim_day)//par->timedInterventions.at(SOCIAL_DISTANCING).at(sim_day)
//             << endl;

        //date,sd,seasonality,vocprev,cinf,closed,rcase,rdeath,Rt
        stringstream ss;
        ss << date->to_string({"yyyy", "mm", "dd"}, "-") << ","
           << community->social_distancing(sim_day)/*par->timedInterventions.at(SOCIAL_DISTANCING).at(sim_day)*/ << ","
           << par->seasonality.at(date->julianDay()-1) << ","
           << (float) community->getNumNewInfections(B_1_1_7)[sim_day]/infections[sim_day] << ","
           << (float) community->getNumNewInfections(B_1_617_2)[sim_day]/infections[sim_day] << ","
           << cAR << ","
           << community->getTimedIntervention(NONESSENTIAL_BUSINESS_CLOSURE, sim_day)<< ","
           << reported_cases*1e4/pop_at_risk << ","
           << rdeaths[sim_day]*1e4/pop_at_risk << ","
           << infections[sim_day]*1e4/pop_at_risk << ","
           << rhosp[sim_day]*1e4/pop_at_risk;
        ledger->plot_log_buffer.push_back(ss.str());

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
//    const vector<size_t> severe             = community->getNumNewlySevere();
//    const vector<size_t> critical           = community->getNumNewlyCritical();
//    const vector<size_t> dead               = community->getNumNewlyDead();
//
    const double cinf  = accumulate(infections.begin(), infections.end(), 0.0);
    const double csymp = accumulate(symptomatic.begin(), symptomatic.end(), 0.0);
//    const double csev  = accumulate(symptomatic.begin(), symptomatic.end(), 0.0);
//    const double ccrit = accumulate(symptomatic.begin(), symptomatic.end(), 0.0);
//    const double cdead = accumulate(symptomatic.begin(), symptomatic.end(), 0.0);
    cerr << "symptomatic infections, total infections, asymptomatic fraction: " << csymp << ", " << cinf << ", " << 1.0 - (csymp/cinf) << endl;
    cerr << "icu deaths, total deaths, ratio: " << cdeath_icu << ", " << cdeath2 << ", " << cdeath_icu/cdeath2 << endl;

//  write_daily_buffer(plot_log_buffer, process_id, "plot_log.csv");
    //return epi_sizes;
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
