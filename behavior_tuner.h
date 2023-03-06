#ifndef _BEHAVIOR_TUNER_H
#define _BEHAVIOR_TUNER_H

using namespace covid::standard;
using namespace covid::util;

class SimulationLedger;
class SimulationCache;

extern gsl_rng* RNG;
extern gsl_rng* REPORTING_RNG;
extern gsl_rng* VAX_RNG;

void gen_simvis(vector<string> &plot_log_buffer);

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
        // community     = nullptr;
        cmty_ledger   = nullptr;
        vc            = nullptr;
        date          = nullptr;
        sim_ledger    = nullptr;
        rng           = nullptr;
        reporting_rng = nullptr;
        vax_rng       = nullptr;
    };

    SimulationCache(Community* o_community, SimulationLedger* o_sim_ledger, gsl_rng* o_rng, gsl_rng* o_reporting_rng, gsl_rng* o_vax_rng) {
        // community     = new Community(*o_community);
        cmty_ledger   = new CommunityLedger(*(o_community->get_ledger()));
        for(Location* hosp : o_community->getLocationsByType(HOSPITAL)) { hosp_people[hosp->getID()] = hosp->getPeople(); }

        vc            = o_community->getVac_Campaign()->quick_cache();

        date          = new Date(*(o_community->get_date()));
        sim_ledger    = new SimulationLedger(*o_sim_ledger);
        rng           = gsl_rng_clone(o_rng);
        reporting_rng = gsl_rng_clone(o_reporting_rng);
        vax_rng       = gsl_rng_clone(o_vax_rng);
    };

    ~SimulationCache() {
        // delete community;
        delete cmty_ledger;
        hosp_people.clear();
        delete vc;
        delete date;
        delete sim_ledger;
        gsl_rng_free(rng);
        gsl_rng_free(reporting_rng);
        gsl_rng_free(vax_rng);
    };

    // Community* community;
    CommunityLedger* cmty_ledger;
    map<int, vector<Person*>> hosp_people;
    Vac_Campaign* vc;
    Date* date;
    SimulationLedger* sim_ledger;
    gsl_rng* rng;
    gsl_rng* reporting_rng;
    gsl_rng* vax_rng;
};

// helps to keep track of the range the binary search can still act on
class Range {
public:
    Range() {};
    Range(double _min, double _max) : min(_min), max(_max) {};

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
        // cerr << "DO YOU WANT TO SIMULATE MANUALLY? (y/n) ";
        // cin >> usr_choice;
        // manual_control = (usr_choice == 'y') ? true :
        //                  (usr_choice == 'n') ? false : true;
        manual_control = false;

        if (not manual_control) {
            // cerr << "DO YOU WANT THE AUTO FITTING TO WAIT FOR KEYPRESSES BEFORE CONTINUING? (y/n) ";
            // cin >> usr_choice;
            // slow_auto = (usr_choice == 'y') ? true :
            //             (usr_choice == 'n') ? false : true;
            slow_auto = false;
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
        output_buffer << right << "  day anchor          window  emp data%  sum dist   cur val  new best              search range   new val";
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

    vector<vector<double>> emp_data;

    stringstream output_buffer;     // buffer to store data to output to the screen during auto tuning
};

// parse empirical data file and store in a map with key of sim_day and value of a vector of cases and deaths
// new req: empty values need to be filled with "NA"
vector<vector<double>> parse_emp_data_file(const Parameters* par) {
    vector<vector<string>> emp_data = read_2D_vector_file(par->rCaseDeathFilename, ',');
    vector<vector<double>> recast_emp_data(2, vector<double>(par->runLength, 0.0));

    for (vector<string> &v : emp_data) {
        if (v[0] == "Date") { continue; }
        const int sim_day = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, v[0]);
        const double rcase = (double) stoi(v[1]);
        // at present our split function doesn't return an empty string element if it's the last element that's missing
        const double rdeath = (v.size() <= 4 or v[4] == "" or v[4] == "NA") ? 0 : (double) stoi(v[4]);

        if (sim_day >= 0 and sim_day < (int) par->runLength) {
            recast_emp_data[0][sim_day] = rcase;
            recast_emp_data[1][sim_day] = rdeath;
        }
    }
    return recast_emp_data;
}

// helper function to handle score the fit in a tuning window
// NOTE: USING DEATHS BY DAY OF DEATH NOT DAY OF REPORTING
double score_fit(const Parameters* par, const Community* community, const size_t sim_day, vector<size_t> sim_rdata_szt, const vector<vector<double>> &emp_data_map) {
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

    // cuts reported cases or deaths from the parsed empirical data
    vector<double> emp_rdata(sim_day + 1, 0.0); // default to 0 if data is missing, for purposes of calculating cumulative sum
    for (size_t day = 0; day <= sim_day; ++day) {
        // if tuning to case data, store cases (index 0), otherwise store deaths (index 1)
        switch (par->behavior_fitting_data_target) {
            case CASES:  emp_rdata[day] = emp_data_map.at(0).at(day); break;
            case DEATHS: emp_rdata[day] = emp_data_map.at(1).at(day); break;
            case NUM_OF_AUTO_FITTING_DATA_TARGETS: [[fallthrough]];
            default: {
                cerr << "ERROR: no PPB auto fitting data target selected." << endl;
                exit(-1);
            }
        }
    }

    const double sim_p10k = 1e4/community->getNumPeople();
    const double fl_pop   = 21538187;
    const double fl_p10k  = 1e4/fl_pop;

    for (auto& val: emp_rdata) { val *= fl_p10k; }     // normalize emp rcases to per 10000 people

    vector<double> sim_rdata(sim_rdata_szt.size());
    for (size_t i = 0; i < sim_rdata.size(); ++i) { sim_rdata[i] = sim_rdata_szt[i] * sim_p10k; }    // normalize sim rcases to per 10000 people

    // calculate cumulative sum of emp and sim rcases
    vector<double> emp_rdata_cumul(emp_rdata.size());
    vector<double> sim_rdata_cumul(sim_rdata.size());
    partial_sum(emp_rdata.begin(), emp_rdata.end(), emp_rdata_cumul.begin());
    partial_sum(sim_rdata.begin(), sim_rdata.end(), sim_rdata_cumul.begin());

    double distance = 0.0;
    double error = 0.0;
    double abs_distance = 0.0;
    double abs_error = 0.0;

    // for each day from the beginning of the sim to now, calculate the error (sim rcases - emp rcases) and the distance (error * weight)
    for (size_t d = 0; d <= sim_day; ++d) {
        //string date = Date::to_ymd(d, par);
        const double daily_crdata_error    = sim_rdata_cumul[d] - emp_rdata_cumul[d];
        const double daily_crdata_distance = daily_crdata_error * fit_weights[d];
        error    += daily_crdata_error;
        distance += daily_crdata_distance;

        abs_error    += abs(daily_crdata_error);
        abs_distance += abs(daily_crdata_distance);
    }

    // handles normalizing distance based on avg rep case incidence (offset by given number of days)
    const size_t fit_norm_offset = par->tuning_window * par->num_preview_windows;   // offset (in days) for when to calculate avg rep case incidence
    size_t num_days_to_avg = 3; // might consider 7 (one week)
    size_t epi_3day_window_start = (int) sim_day - fit_norm_offset - (num_days_to_avg - 1)/2;   // day to start averaging at
    vector<size_t> rdata_to_avg(emp_rdata.begin() + epi_3day_window_start, emp_rdata.begin() + epi_3day_window_start + num_days_to_avg);
    const double mean_rdata  = mean(rdata_to_avg);

    vector<double> single_rdata_p10k = {0, 0, (1 * fl_p10k)};   // for use if there are no rep cases to average (needed to prevent division by 0 later on)

    double avg_3day_rdata_w_offset = mean_rdata == 0 ? mean(single_rdata_p10k) : mean_rdata;

    const double normed_distance = distance/avg_3day_rdata_w_offset;
    return normed_distance;
}

// helper function to conduct modified binary search
double bin_search_anchor(BehaviorAutoTuner* tuner, double distance) {
    // grab the current search range and anchor value from the tuner
    Range* range = tuner->bin_search_range;
    const double cur_val = tuner->cur_anchor_val;

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
BehaviorAutoTuner* initialize_behavior_autotuning(const Parameters* par) {
    BehaviorAutoTuner* tuner = new BehaviorAutoTuner();
    tuner->init_defaults();
    tuner->user_sim_setup();

    vector<vector<double>> parsed_data = parse_emp_data_file(par);
    vector<double> case_ma  = calc_centered_avg(parsed_data[0], 7);
    vector<double> death_ma = calc_centered_avg(parsed_data[1], 7);
    tuner->emp_data = {case_ma, death_ma};

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
    sim_cache = new SimulationCache(community, ledger, RNG, REPORTING_RNG, VAX_RNG);
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

        double val_next;
        if (tuner->manual_control) {
            cerr << "FIT IS GOOD" << endl << "Enter next soc_contact step: ";
            cin >> val_next;
        } else {
            val_next = tuner->cur_anchor_val;   // always start the next window with the val that was selected
            tuner->reset_bin_search_range();
            if (tuner->slow_auto) { cin.ignore(); }

            tuner->output_buffer << right << setw(10) << tuner->cur_anchor_val;
            cerr << right << tuner->output_buffer.str() << endl;
            tuner->clear_output_buffer();
//            cerr << endl << endl;
//            tuner->print_header();
        }

        // need to overwrite the previously changing anchor with the best performing value
        if (tuner->tuning_window_ct == 1) { // tuning the first window
            double val1, val2;
            if (tuner->manual_control) {
                // TODO: fix manual mode
                exit(-1);
            } else {
                val1 = tuner->cur_anchor_val;
                val2 = val1;
            }

            social_distancing_anchors[0] = {Date::to_ymd(0, par), val1};
            social_distancing_anchors[1] = {Date::to_ymd(par->tuning_window - 1, par), val2};
        } else { // past the first tuning window
            double val;
            if (tuner->manual_control) {
                // TODO: fix manual mode
                exit(-1);
            } else {
                val = tuner->cur_anchor_val;
            }

            const size_t anchor_to_update_day = (par->tuning_window * tuner->tuning_window_ct) - 1;
            social_distancing_anchors.back() = {Date::to_ymd(anchor_to_update_day, par), val};
        }

        const size_t next_anchor_day = (par->tuning_window * (tuner->tuning_window_ct + 1)) - 1;
        social_distancing_anchors.emplace_back(Date::to_ymd(next_anchor_day, par), val_next);
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

            const size_t anchor_to_update_day = (par->tuning_window * tuner->tuning_window_ct) - 1;
            social_distancing_anchors.back() = {Date::to_ymd(anchor_to_update_day, par), val};
        }
        //  this is what "waits" for the 'enter' keypress
        if (not tuner->manual_control and tuner->slow_auto) { cin.ignore(); }

        // output tuning data to the screen
//        cerr << right << tuner->output_buffer.str() << endl;
        tuner->clear_output_buffer();
     }
}


bool restore_from_cache(Community* &community, Date* &date, SimulationCache* sim_cache, SimulationLedger* &ledger, vector<TimeSeriesAnchorPoint> social_distancing_anchors) {
    if (ledger) { delete ledger; }
    ledger = new SimulationLedger(*(sim_cache->sim_ledger));

    // if (community) { delete community; }
    // community = new Community(*(sim_cache->community));
    community->load_from_cache(sim_cache->cmty_ledger, sim_cache->date, sim_cache->hosp_people, sim_cache->vc);
    date      = community->get_date();
    community->setSocialDistancingTimedIntervention(social_distancing_anchors);

    if (sim_cache->rng)           { gsl_rng_memcpy(RNG, sim_cache->rng); }
    if (sim_cache->reporting_rng) { gsl_rng_memcpy(REPORTING_RNG, sim_cache->reporting_rng); }
    if (sim_cache->vax_rng)       { gsl_rng_memcpy(VAX_RNG, sim_cache->vax_rng); }

    if (community->getVac_Campaign()) { community->getVac_Campaign()->set_rng(VAX_RNG); }
    return true;
}


// main helper function to control the behavior auto tuning system
void behavior_autotuning(const Parameters* par, Community* &community, Date* &date, SimulationLedger* &ledger, BehaviorAutoTuner* tuner, SimulationCache* &sim_cache, vector<TimeSeriesAnchorPoint> &social_distancing_anchors, bool& restore_occurred) {
    const size_t day = date->day();
    const size_t recaching_day = (tuner->tuning_window_ct * par->tuning_window) - 1;
    size_t behavior_processing_day;
    switch (par->behavior_fitting_data_target) {
        case CASES:  behavior_processing_day = ((tuner->tuning_window_ct + par->num_preview_windows) * par->tuning_window) - 1; break;
        case DEATHS: behavior_processing_day = (((tuner->tuning_window_ct + par->num_preview_windows) * par->tuning_window) - 1) + par->death_tuning_offset; break;
        case NUM_OF_AUTO_FITTING_DATA_TARGETS: [[fallthrough]];
        default: {
            cerr << "ERROR: no PPB auto fitting data target selected." << endl;
            exit(-1);
        }
    }


    if (tuner->recache and (day == recaching_day)) {
        // we are at the end of tuning window that was deemed "good"
        overwrite_sim_cache(sim_cache, community, ledger, tuner);   // TODO: change to quick cache overwrite
    } else if (day == behavior_processing_day) {
        // calculates the proportion of days in the tuning window that has empirical data
        size_t window_start_sim_day = (day + 1) - ((par->num_preview_windows + 1) * par->tuning_window);
        //size_t num_days_w_emp_data = 0;
        //for (size_t d = window_start_sim_day; d <= day; ++d) { if (tuner->emp_data.count(d)) { ++num_days_w_emp_data; } }
        double prop_days_w_emp_data = 1.0;//(((double) num_days_w_emp_data) / (day - window_start_sim_day + 1));

        // begin adding tuning data to print to screen
        tuner->output_buffer << right
                             << setw(5) << day
                             << setw(7) << (tuner->tuning_window_ct * par->tuning_window) - 1
                             << setw(3) << "[" << setw(5) << window_start_sim_day << ", " << setw(5) << day << "]"
                             << setw(10) << prop_days_w_emp_data * 100 << "%";

        // calculate normalized distance for this tuning window
        vector<size_t> model_tuning_data;
        switch (par->behavior_fitting_data_target) {
            case CASES: model_tuning_data = community->getNumDetectedCasesReport(); break;
            case DEATHS: model_tuning_data = community->getNumDetectedDeathsReport(); break;
            case NUM_OF_AUTO_FITTING_DATA_TARGETS: [[fallthrough]];
            default: {
                cerr << "ERROR: no PPB auto fitting data target selected." << endl;
                exit(-1);
            }
        }
        double fit_distance = score_fit(par, community, day, model_tuning_data, tuner->emp_data);
        //cerr << "current anchor, score: " << tuner->cur_anchor_val << ", " << fit_distance << endl;
        // keep track of the anchor val with the smallest distance for this window
        tuner->update_best_distance(fit_distance);

//        gen_simvis(ledger->plot_log_buffer);

        int fit_is_good = 0;
        if (tuner->manual_control) {
            cerr << "IS THE FIT GOOD? ";
            cin >> fit_is_good;
        } else {
                                                                // is the abs(distance) less than the threshold after adjusting for the presence of emp data
            fit_is_good = abs(fit_distance) < (tuner->fit_threshold * prop_days_w_emp_data)
                                                                // is the distance pos and the cur val is the max of the search range (can't search more)
                          or (fit_distance > 0 and tuner->cur_anchor_val == tuner->bin_search_range_max())
                                                                // is the distance neg and the cur val is the min of the search range (can't search more)
                          or (fit_distance < 0 and tuner->cur_anchor_val == tuner->bin_search_range_min())
                                                                // are we out of data to fit against
                          or (tuner->tuning_window_ct > 1 and prop_days_w_emp_data == 0)
                                                                // are the search range min and max equal (can't search more)
                          or (tuner->bin_search_range_max() == tuner->bin_search_range_min())
                                                                // is this before people started reacting 
                          or *date <= "2020-03-01";

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
        restore_occurred = restore_from_cache(community, date, sim_cache, ledger, social_distancing_anchors);   // TODO: change to quick revert
    }
}


// if using previous tuned values, parse that file and save the behavior anchors
// anchor file expected to be csv with columns date, anchor_val
void init_behavioral_vals_from_file(const Parameters* par, Community* community) {
    vector< vector<string> > tuning_data = read_2D_vector_file(par->behaviorInputFilename, ',');
    vector<TimeSeriesAnchorPoint> tuned_anchors;

    for (vector<string> &v : tuning_data) {
        assert(v.size() == 2);
        if (v[0] == "date") { continue; }
        tuned_anchors.emplace_back(v[0], stod(v[1]));
    }

    community->setSocialDistancingTimedIntervention(tuned_anchors);
}

void write_anchors_to_file(const Parameters* par, vector<TimeSeriesAnchorPoint> anchors) {
    ofstream ofs(par->behaviorOutputFilename);
    ofs << "date,anchor_val" << endl;
    for (const TimeSeriesAnchorPoint &tsap : anchors) {
        ofs << tsap.date << "," << setprecision(20) << tsap.value << setprecision(6) << endl;
    }
    ofs.close();
}

#endif
