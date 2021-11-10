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

// TODO(alex): should be moved to Parameters.h eventually
enum FittingDataSource { EMP, SIM, NUM_OF_FITTING_DATA_SOURCES };
// TODO(alex): move somewhere more proper
gsl_rng* FITTING_RNG = gsl_rng_alloc(gsl_rng_mt19937);

/*
HEURISTIC PLANNING:

SCENARIOS (assume that the window is already deemed to have a bad fit):
- all days have positive error --> increase param nearest day with max abs(error)
- all days have negative error --> decrease param nearest day with max abs(error)
- error changes sign mid-window
    - positive error > negative error --> increase param nearest day with max error
    - positive error < negative error --> decrease param nearest day with min error

*/

/*size_t closest_window_param(const signed int day_to_check, const signed int window_start, const signed int window_mid, const signed int window_end) {
//    map<int, size_t> day_to_check_diffs;
//    vector<int> window_diffs;
//
//    day_to_check_diffs[abs(day_to_check - window_start)] = (size_t) window_start;
//    window_diffs.push_back(abs(day_to_check - window_start));
//
//    day_to_check_diffs[abs(day_to_check - window_mid)]   = (size_t) window_mid;
//    window_diffs.push_back(abs(day_to_check - window_mid));
//
//    day_to_check_diffs[abs(day_to_check - window_end)]   = (size_t) window_end;
//    window_diffs.push_back(abs(day_to_check - window_end));
//
//    return day_to_check_diffs[min_element(window_diffs)];

    if (day_to_check > window_mid) {
        return (size_t) window_mid;
    } else {
        return (size_t) window_start;
    }
}*/

/*pair<size_t, double> social_contact_param_to_adj(const Parameters* par, const size_t sim_day, map<size_t, TimeSeriesAnchorPoint> &social_contact_map,
                                                 const map<size_t, double> &window_fit_map) {
    size_t window_end = 0;
    for (size_t adj = 0; adj < par->fitting_window; ++adj) {
        if ((sim_day + adj + 1) % par->fitting_window == 0) {
            window_end = (signed int) (sim_day + adj);
        }
    }
    const size_t window_start = window_end - par->fitting_window + 1;
    const size_t window_mid   = window_start + (par->fitting_window/2);

    // DO HEURISTIC SEARCH HERE
    unsigned int num_pos_error = 0, num_neg_error = 0;
    size_t max_abs_error_day = 0, max_error_day = 0, min_error_day = 0;
    double max_abs_error = 0, max_error = INT_MIN, min_error = INT_MAX;
    double sum_error = 0.0;
    for (const auto& [day, error] : window_fit_map) {
        (error > 0) ? num_pos_error++ : num_neg_error++;
        sum_error += error;

        if (abs(error) > max_abs_error) {
            max_abs_error = abs(error);
            max_abs_error_day   = (signed int) day;
        }

        if (error > max_error) {
            max_error = error;
            max_error_day   = (signed int) day;
        }

        if (error < min_error) {
            min_error = error;
            min_error_day   = (signed int) day;
        }
    }

    size_t day_to_adj = window_end;

    if ((not num_pos_error == 0) and (num_neg_error == 0)) {
        // all days have pos error --> increase param nearest day with max abs(error)
        day_to_adj = closest_window_param((signed int) max_abs_error_day, (signed int) window_start, (signed int) window_mid, (signed int) window_end);
    } else if ((num_pos_error == 0) and (not num_neg_error == 0)) {
        // all days have neg error --> decrease param nearest day with max abs(error)
        day_to_adj = closest_window_param((signed int) max_abs_error_day, (signed int) window_start, (signed int) window_mid, (signed int) window_end);
    } else if ((not num_pos_error == 0) and (not num_neg_error == 0)) {
        // sign changes mid-window
        if (num_pos_error > num_neg_error) {
            // positive error > negative error --> increase param nearest day with max error
            day_to_adj = closest_window_param((signed int) max_error_day, (signed int) window_start, (signed int) window_mid, (signed int) window_end);
        } else if (num_pos_error < num_neg_error) {
            // positive error < negative error --> decrease param nearest day with min error
            day_to_adj = closest_window_param((signed int) min_error_day, (signed int) window_start, (signed int) window_mid, (signed int) window_end);
        }
    }

    if (social_contact_map.count(day_to_adj)) {
        // param exists at the select day
  cerr << "PARAM AT " << day_to_adj << "(" << social_contact_map[day_to_adj].value << ")" << " >1? " << (social_contact_map[day_to_adj].value >= 1) << endl;
        if ((social_contact_map[day_to_adj].value <= 0) or (social_contact_map[day_to_adj].value >= 1)) {
            // can't be adjusted any further
  cerr << "CANT ADJUST " << day_to_adj << "(" << social_contact_map[day_to_adj].value << ")" << endl;
            while (not (day_to_adj == window_start)) {
                // cycle through the other available parameters before this parameter
                day_to_adj = (day_to_adj == window_end) ? window_mid :
                             (day_to_adj == window_mid) ? window_start : window_start;

                if (social_contact_map.count(day_to_adj)) {
                    // param exists at the selected day
  cerr << "PARAM AT " << day_to_adj << endl;
                    if ((social_contact_map[day_to_adj].value <= 0) or (social_contact_map[day_to_adj].value >= 1)) {
                        // selected day can't be adjusted any further
  cerr << "CANT ADJUST " << day_to_adj << "(" << social_contact_map[day_to_adj].value << ")" << endl;
                        continue;
                    } else {
                        break;
                    }
                }
            }
            if ((day_to_adj == window_start) and (((social_contact_map[day_to_adj].value <= 0) or (social_contact_map[day_to_adj].value >= 1)))) {
                // all params in the window can not be altered furthered
                cerr << "NO WINDOW PARAMETERS TO ADJUST" << endl;
                exit(-222);
            }
        }
    }

cerr << "DAY TO ADJ: " << day_to_adj << endl;
    // DO NEW PARAM CALCULATION HERE
    // logit transform existing value
cerr << "ORIG PARAM: " << social_contact_map[day_to_adj].value << endl;
    double logit_transformed_param = logit(social_contact_map[day_to_adj].value);
cerr << "LOGIT TRANSFORM: " << logit_transformed_param << endl;
    // select adjustment from gaussian distribution (with stdev informed by the fit error)
    double param_adjustment = (sum_error > 0) ? abs(gsl_ran_gaussian(FITTING_RNG, (sum_error/1e1))) : -abs(gsl_ran_gaussian(FITTING_RNG, (sum_error/1e1)));
cerr << "ADJUSTMENT: " << param_adjustment << endl;
    // adjust the parameter and transform back into logistic space
    double new_param = logistic(logit_transformed_param + param_adjustment);
cerr << "NEW PARAM: " << new_param << endl;

    return pair<size_t, double>{day_to_adj, new_param};
}*/

/*pair<bool, pair<size_t, double> > check_window_social_contact_fit(const Parameters* par, map<size_t, int> &emp_data, const vector<size_t> &sim_rcases, const size_t sim_day_today,
                                                   map<size_t, TimeSeriesAnchorPoint> &social_contact_map, const double fitting_threshold) {
    const size_t emp_start = emp_data.begin()->first;
    const size_t emp_end   = emp_data.end()->first;
    const size_t window_start = (sim_day_today + 1) - par->fitting_window;
    // if the emp data does not contain vlaues for the sim window, continue the simulation
    if ( ((window_start < emp_start) and (sim_day_today < emp_start)) or
         ((window_start > emp_end)   and (sim_day_today > emp_end)) ) {
         cerr << "NO EMP DATA TO FIT TO" << endl; 
         return pair<bool, pair<size_t, double> >{true, {par->runLength, 0.0}};
    } else { // check fit here: return true if fit is acceptable and false to re-simulate
        const double fl_per_10k = 1e4/20609673;//21538187; // TODO(alex): i dont like having these numbers here, but dont know how to solve it yet; will come back to this later
        const double sim_per_10k = 1e4/375474;

        const size_t overlap_start = max(emp_start, window_start);
        const size_t overlap_end   = min(emp_end, sim_day_today);

        // work with only the data for up to and including the simulated window of time
        vector< vector<double> > relevant_data(NUM_OF_FITTING_DATA_SOURCES, vector<double>(overlap_end + 1, 0.0));
        for (size_t day = 0; day <= overlap_end; ++day) {
            relevant_data[EMP].at(day) = (day >= emp_start) ? emp_data.at(day) * fl_per_10k : 0;
            relevant_data[SIM].at(day) = sim_rcases.at(day) * sim_per_10k;
        }

        // calculate cum sum of rcases for emp and sim
        vector< vector<double> > cumsum_data(NUM_OF_FITTING_DATA_SOURCES, vector<double>(overlap_end + 1, 0.0));
        for (size_t i = 0; i < NUM_OF_FITTING_DATA_SOURCES; ++i) {
            double cumsum = 0.0;
            for (size_t day = 0; day <= overlap_end; ++day) {
                cumsum += relevant_data[i].at(day);
                cumsum_data[i].at(day) = cumsum;
            }
        }

        // calculate abs and rel error for each day in the overlap window
        cerr << "DAY\t\tSIM RCs\t\tEMP RCs\t\tABS ERROR\t5\% EMP" << endl;
        vector<double> abs_error, rel_error;
        map<size_t, double> window_fit_map;
        for (size_t i = overlap_start; i <= overlap_end; ++i) {
            const double daily_abs_error = cumsum_data[SIM].at(i) - cumsum_data[EMP].at(i);
            const double daily_rel_error = daily_abs_error / cumsum_data[EMP].at(i);
            const double fitting_boundary = cumsum_data[EMP].at(i) * 0.05;
            if (daily_rel_error == -1.0) { continue; }
            abs_error.push_back(daily_abs_error);
            window_fit_map[i] = daily_abs_error;
            rel_error.push_back(daily_rel_error);
            cerr << i << "\t\t" << cumsum_data[SIM].at(i) << "\t\t" << cumsum_data[EMP].at(i) << "\t\t" << daily_abs_error << "\t\t" << fitting_boundary << endl;
        }

        if (abs_error.size() and rel_error.size()) {
            const double total_window_abs_error = sum(abs_error);
            const double avg_window_abs_error   = mean(abs_error);
            const double avg_window_rel_error   = mean(rel_error);
            cerr << "\nSUM OF ABS ERROR: " << total_window_abs_error << "\nAVG OF ABS ERROR: " << avg_window_abs_error << "\nAVG OF REL ERROR: " << avg_window_rel_error << "\n" << endl;

            if (abs(total_window_abs_error) <= (fitting_threshold * par->fitting_window)) {
                cerr << "ACCEPTABLE FIT" << endl;
                return pair<bool, pair<size_t, double> >{true, {par->runLength, 0.0}};
            } else {
                pair<size_t, double> fit_adjustment = social_contact_param_to_adj(par, sim_day_today, social_contact_map, window_fit_map);
                cerr << "UNACCEPTABLE FIT" << endl;
                return pair<bool, pair<size_t, double> >{false, fit_adjustment};
            }
 
        } else {
            // we are here if there were no days to make meaningful comparisons between
            cerr << "NOTHING TO COMPARE" << endl;
            return pair<bool, pair<size_t, double> >{true, {par->runLength, 0.0}};
        }
    }
}*/

// pair<bool, double> check_daily_social_contact_fit(const Parameters* par, map<size_t, int> &emp_data, const vector<size_t> &sim_rcases, const size_t sim_day_today,
//                                                      const double fitting_threshold) {
//     const size_t emp_start = emp_data.begin()->first;
//     const size_t emp_end   = emp_data.end()->first;
// 
//     if ((sim_day_today < emp_start) or(sim_day_today > emp_end)) { // if the emp data does not contain vlaues for the sim window, continue the simulation
//         //cerr << "NO EMP DATA TO FIT TO" << endl;
//         return pair<bool, double>{true, 0.0};
//     } else { // check fit here: return true if fit is acceptable and false to re-simulate
//         const double fl_per_10k = 1e4/20609673;//21538187; // TODO(alex): i dont like having these numbers here, but dont know how to solve it yet; will come back to this later
//         const double sim_per_10k = 1e4/375474;
// 
//         // work with only the data for up to and including the simulated window of time
//         vector< vector<double> > relevant_data(NUM_OF_FITTING_DATA_SOURCES, vector<double>(sim_day_today + 1, 0.0));
//         for (size_t day = 0; day <= sim_day_today; ++day) {
//             relevant_data[EMP].at(day) = ((day >= emp_start) and (day <= emp_end)) ? emp_data.at(day) * fl_per_10k : 0;
//             relevant_data[SIM].at(day) = sim_rcases.at(day) * sim_per_10k;
//         }
// 
//         // calculate cum sum of rcases for emp and sim
//         vector< vector<double> > cumsum_data(NUM_OF_FITTING_DATA_SOURCES, vector<double>(sim_day_today + 1, 0.0));
//         for (size_t i = 0; i < NUM_OF_FITTING_DATA_SOURCES; ++i) {
//             double cumsum = 0.0;
//             for (size_t day = 0; day <= sim_day_today; ++day) {
//                 cumsum += relevant_data[i].at(day);
//                 cumsum_data[i].at(day) = cumsum;
//             }
//         }
// 
//         // calculate abs and rel error for each day in the overlap window
//         const double today_abs_error = cumsum_data[SIM].at(sim_day_today) - cumsum_data[EMP].at(sim_day_today);
//         const double today_rel_error = today_abs_error / cumsum_data[EMP].at(sim_day_today);
// 
//         if (abs(today_abs_error) > 0 and not (today_rel_error == -1)) {
// cerr << "DAY\t\tSIM RCs\t\tEMP RCs\t\tABS ERROR\tREL ERROR" << endl;
// cerr << sim_day_today << "\t\t" << cumsum_data[SIM].at(sim_day_today) << "\t\t" << cumsum_data[EMP].at(sim_day_today) << "\t\t" << today_abs_error << "\t\t" << today_rel_error << endl;
//             const double fitting_boundary = cumsum_data[EMP].at(sim_day_today) * fitting_threshold;
// cerr << "\nTODAY ABS ERROR: " << today_abs_error << "\nFITTING BOUNDARY: " << fitting_boundary << endl;
//             if (abs(today_abs_error) > fitting_boundary) {
//                 // bad fit
// cerr << "BAD FIT " << today_abs_error << endl << endl;
//                 return pair<bool, double>{false, today_abs_error};
//             } else {
//                 // good fit
// cerr << "GOOD FIT" << endl << endl;
//                 return pair<bool, double>{true, 0.0};
//             }
//         } else {
//             // we are here if there were no days to make meaningful comparisons between
//             return pair<bool, double>{true, 0.0};
//         }
//     }
// }

/*void alter_social_contact_params(const Parameters* par, const pair<size_t, double> adjustment, map<size_t, TimeSeriesAnchorPoint> &social_contact_map) {
    const size_t day_to_adj               = get<size_t>(adjustment);
    const double social_contact_param_adj = get<double>(adjustment);

    Date* tmp_date = new Date(par);
    for(; (size_t) tmp_date->day() <= day_to_adj; tmp_date->increment()) {
        if ((size_t) tmp_date->day() == day_to_adj) {
            if (social_contact_map.count(day_to_adj)) {
                social_contact_map[day_to_adj].value = social_contact_param_adj;
            } else {
                TimeSeriesAnchorPoint tsap = {tmp_date->to_ymd(), social_contact_param_adj};
                social_contact_map[day_to_adj] = tsap;
            }

            if (social_contact_map[day_to_adj].value > 0.99)      { social_contact_map[day_to_adj].value = 1; }
            else if (social_contact_map[day_to_adj].value < 0.01) { social_contact_map[day_to_adj].value = 0; }
        }
    }
    delete tmp_date;
}*/

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
    int retval = system("Rscript fitvis.R");
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

    double min;
    double max;
};

double bin_search_anchor(Range &range, double cur_val, double distance) {
    const double MIN_ADJ = 0.01;
    double new_val;
cerr << "CUR BIN SEARCH RANGE: [" << range.min << ", " << range.max << "]" << endl;
cerr << "CUR VAL, fit distance:              " << cur_val << ", " << distance << endl;
    assert(cur_val <= range.max);
    assert(cur_val >= range.min);
    assert(range.max >= range.min);

    range = distance > 0 ? Range(cur_val, range.max) : Range(range.min, cur_val);
    const double adj = (range.max - range.min) / 2.0;
    if (adj < MIN_ADJ) {
        if (distance > 0) {
            new_val = range.max;
            range.min = range.max;
        } else {
            new_val = range.min;
            range.max = range.min;
        }
    } else {
        new_val = range.min + adj;
    }
cerr << "ADJ, NEW VAL:         " << adj << ", " << new_val << endl;
cerr << "NEW BIN SEARCH RANGE: [" << range.min << ", " << range.max << "]" << endl;
    return new_val;
}

vector<string> simulate_epidemic(const Parameters* par, Community* &community, const string /*process_id*/, const vector<string> mutant_intro_dates) {//,
                                 //map<size_t, TimeSeriesAnchorPoint> &social_contact_map) {
    SimulationLedger* ledger    = new SimulationLedger();
    SimulationCache* sim_cache  = nullptr;

    //bool window_fit_is_good = true;
    //vector<int> epi_sizes;
    //vector<string> daily_output_buffer;

    //map<string, vector<int> > periodic_incidence = construct_tally();
    ledger->periodic_incidence = construct_tally();
    //vector<int> periodic_prevalence(NUM_OF_PREVALENCE_REPORTING_TYPES, 0);
    ledger->periodic_prevalence = vector<int>(NUM_OF_PREVALENCE_REPORTING_TYPES, 0);
    vector<double> trailing_averages(par->runLength);
    const double pop_at_risk = min(community->getNumPeople(), par->numSurveilledPeople);
    vector<TimeSeriesAnchorPoint> social_distancing_anchors; /*= {
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
    };*/
    //community->setSocialDistancingTimedIntervention(social_distancing_anchors);

    //vector<string> plot_log_buffer = {"date,sd,seasonality,vocprev1,vocprev2,cinf,closed,rcase,rdeath,inf,rhosp,Rt"};
    ledger->plot_log_buffer = {"date,sd,seasonality,vocprev1,vocprev2,cinf,closed,rcase,rdeath,inf,rhosp,Rt"};

    Date* date = community->get_date();
    ledger->strains = {50.0, 0.0, 0.0}; // initially all WILDTYPE
    assert(ledger->strains.size() == NUM_OF_STRAIN_TYPES);

    map<size_t, vector<int>> emp_data = parse_emp_data_file(par, "rcasedeath-florida.csv");

    size_t fitting_window_ct = 1;

    ofstream window_fit_file;
    window_fit_file.open("window_fit_dates.out", ofstream::out | ofstream::trunc);
    window_fit_file << "date,val" << endl;

    ofstream fit_err_file;
    fit_err_file.open("fit_error.out", ofstream::out | ofstream::trunc);
    fit_err_file << "date,abs_err,adj_err,norm_abs_err,norm_adj_err" << endl;


    Range bin_search_range(0.0, 1.0);
    double cur_anchor_val  = 0.0;

    cerr << "DO YOU WANT TO SIMULATE MANUALLY? (y/n) ";
    char usr_choice;
    cin >> usr_choice;
    bool manual_control = (usr_choice == 'y') ? true :
                          (usr_choice == 'n') ? false : true;

    bool slow_auto = true;
    if (not manual_control) {
        cerr << "DO YOU WANT THE AUTO FITTING TO WAIT FOR KEYPRESSES BEFORE CONTINUING? (y/n) ";
        cin >> usr_choice;
        slow_auto = (usr_choice == 'y') ? true :
                    (usr_choice == 'n') ? false : true;
    }

    if (par->auto_fitting) {
        // inital sim cache
        sim_cache = new SimulationCache(community, ledger, RNG, REPORTING_RNG);

        community->_clearSocialDistancingTimedIntervention();
        cerr << "ENTER SOC_CONTACT PARAMS FOR FIRST AND SECOND ANCHOR POINTS: ";

        double val1, val2;
        if (manual_control) {
            cin >> val1 >> val2;
        } else {
            val1 = 0.0;//0.25;
            val2 = val1;
            cur_anchor_val  = val1;
        }
        social_distancing_anchors.emplace_back(Date::to_ymd(0, par), val1);
        social_distancing_anchors.emplace_back(Date::to_ymd(par->fitting_window - 1, par), val2);
        community->setSocialDistancingTimedIntervention(social_distancing_anchors);

        for (auto anchor: social_distancing_anchors) { window_fit_file << anchor.date << "," << anchor.value << endl; }
    }

    bool recache = false;

    const double fit_threshold = 20.0;//25.0;//50.0;

    for (; date->day() < (signed) par->runLength; date->increment()) {
        community->tick();

        if (par->auto_fitting) {
            const size_t day = date->day();
            if (recache and (day + 1) == fitting_window_ct * par->fitting_window) {
cerr << "RECACHING ON DAY " << day + 1 << endl;
                if (sim_cache) { delete sim_cache; } 
                sim_cache = new SimulationCache(community, ledger, RNG, REPORTING_RNG);
                fitting_window_ct++;
                recache = false;
            } else if ((day + 1) == (fitting_window_ct + par->num_preview_windows) * par->fitting_window) {
size_t window_start_sim_day = (day + 1) - ((par->num_preview_windows + 1) * par->fitting_window);
cerr << "FITTING WINDOW: " << window_start_sim_day << " TO " << window_start_sim_day + par->fitting_window - 1 << endl;
cerr << "PREVIEW WINDOW: " << window_start_sim_day + par->fitting_window  << " TO " << day << endl;
                vector<double> fit_error = fitting_error(par, community, day, community->getNumDetectedCasesReport(), emp_data);


                gen_simvis(ledger->plot_log_buffer);
//                for (auto anchor: social_distancing_anchors) { cerr << anchor.date << ": " << anchor.value << endl; }
                int tmp_cnt = 0;
                for (auto anchor: social_distancing_anchors) {
                    if (tmp_cnt == 0) { cerr << anchor.value << " "; tmp_cnt++; }
                    else { cerr << anchor.value << " 1 "; }
                }
                cerr << endl;
                // for (size_t d = 0; d <= day; ++d) { cerr << d << ": " << community->social_distancing(d) << endl; }
                cerr << "IS THE FIT GOOD? ";
                int fit_is_good = 0;

                if (manual_control) {
                    cin >> fit_is_good;
                } else {
                    cerr << endl;
                    fit_is_good = abs(fit_error[1]) < fit_threshold
                                  or (fit_error[1] > 0 and cur_anchor_val == bin_search_range.max)
                                  or (fit_error[1] < 0 and cur_anchor_val == bin_search_range.min)
                                  or (bin_search_range.max == bin_search_range.min);
                }

                if(fit_is_good) {
           cerr << "FIT IS GOOD" << endl;
                    recache = true;
                    cerr << "Enter next soc_contact step: ";

                    double val1;
                    if (manual_control) {
                        cin >> val1;
                    } else {
                        cerr << cur_anchor_val << endl;
                        val1 = cur_anchor_val;
                        bin_search_range = {0.0, 1.0};
                        if (slow_auto) { cin.ignore(); }
                    }
                    social_distancing_anchors.emplace_back(Date::to_ymd((par->fitting_window * (fitting_window_ct + 1)) - 1, par), val1);

                    for (auto anchor : social_distancing_anchors) { window_fit_file << anchor.date << "," << anchor.value << endl; }
                    fit_err_file << Date::to_ymd(window_start_sim_day + par->fitting_window - 1, par) << "," << fit_error[0] << "," << fit_error[1] << "," 
                                 << -1 << "," << -1 << endl;
                } else {
          cerr << "FIT IS NOT GOOD" << endl;
                    if (not manual_control) { cur_anchor_val = bin_search_anchor(bin_search_range, cur_anchor_val, fit_error[1]); }
                    if (fitting_window_ct == 1) {
                        cerr << "Enter soc_contact params for first and second anchor points: ";

                        double val1, val2;
                        if (manual_control) {
                            cin >> val1 >> val2;
                        } else {
                            cerr << cur_anchor_val << endl;
                            val1 = cur_anchor_val;
                            val2 = val1;
                        }
                        social_distancing_anchors[0] = {Date::to_ymd(0, par), val1};
                        social_distancing_anchors[1] = {Date::to_ymd(par->fitting_window - 1, par), val2};

                        for (auto anchor: social_distancing_anchors) { window_fit_file << anchor.date << "," << anchor.value << endl; }
                    } else {
                        cerr << "Enter new value for this window's soc_contact param: ";
                        double val;
                        if (manual_control) {
                            cin >> val;
                        } else {
                            cerr << cur_anchor_val << endl;
                            val = cur_anchor_val;
                        }
                        social_distancing_anchors.back() = {Date::to_ymd((par->fitting_window * fitting_window_ct) - 1, par), val};

                        for (auto anchor: social_distancing_anchors) { window_fit_file << anchor.date << "," << anchor.value << endl; }
                    }
                    if (not manual_control and slow_auto) { cin.ignore(); }
                }

                if (ledger) { delete ledger; }
                ledger = new SimulationLedger(*(sim_cache->sim_ledger));

                if (community) { delete community; }
                community = new Community(*(sim_cache->community));
                date      = community->get_date();
                community->setSocialDistancingTimedIntervention(social_distancing_anchors);

                if (sim_cache->rng)           { gsl_rng_memcpy(RNG, sim_cache->rng); }
                if (sim_cache->reporting_rng) { gsl_rng_memcpy(REPORTING_RNG, sim_cache->reporting_rng); }
            }
        }

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
        const vector<size_t> rdeaths            = community->getNumDetectedDeaths();
        //const vector<size_t> rdeaths            = community->getNumNewlyDead();


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
    if (window_fit_file.is_open()) { window_fit_file.close(); }
    if (fit_err_file.is_open())    { fit_err_file.close(); }

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
    if (FITTING_RNG) { gsl_rng_free(FITTING_RNG); }

    vector<string> plot_log_buffer = ledger->plot_log_buffer;
    if (ledger) { delete ledger; }

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
