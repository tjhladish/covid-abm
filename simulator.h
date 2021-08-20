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

/*void update_vaccinations(const Parameters* par, Community* community, const Date* date) {
    const int doseInterval = par->vaccineDoseInterval;
    assert(doseInterval > 0); // only non-zero positive ints are sensible
    //const int boostInterval = par->vaccineBoostingInterval;
    for (CatchupVaccinationEvent cve: par->catchupVaccinationEvents) {
        // Normal, initial vaccination -- boosting, multiple doses handled in Community::tick()
        if (date->day() >= (signed) cve.campaignStart and date->day() < (signed) (cve.campaignStart + cve.campaignDuration)) {
            if (par->abcVerbose) cerr << "vaccinating " << cve.coverage*100 << "% of age " << cve.age << " over " << cve.campaignDuration << " days starting on day " << cve.campaignStart << endl;
            community->vaccinate(cve);
        }
    }
}*/


int seed_epidemic(const Parameters* par, Community* community, const Date* date, vector<StrainType> strains) {
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
            if (community->infect(transmit_to_id, choice(RNG, strains))) {
                introduced_infection_ct++;
            }
        }
    }
    return introduced_infection_ct;
}


void advance_simulator(const Parameters* par, Community* community, Date* date, const string process_id, map<string, vector<int> > &periodic_incidence, vector<int> &periodic_prevalence, vector<int> &epi_sizes) {
    community->tick();

    vector<StrainType> strains = {WILDTYPE};
    seed_epidemic(par, community, date, strains);

    for (Person* p: community->getPeople()) {
        const int now = date->day();

        if (p->isInfected(now) or p->isSymptomatic(now) or p->isDead(now)) {
            const Infection* infec = p->getInfection();
            bool intro  = not infec->isLocallyAcquired();

            // Incidence
            if (infec->getInfectedTime() == now) {
                periodic_incidence["daily"][INTRO_INF]      += intro;
                periodic_incidence["daily"][TOTAL_INF]      += 1;
            }

            if (infec->getSymptomTime() == now) {
                periodic_incidence["daily"][INTRO_CASE]     += intro;
                periodic_incidence["daily"][TOTAL_CASE]     += 1;
            }

            if (infec->getSevereTime() == now) {
                periodic_incidence["daily"][INTRO_SEVERE]   += intro;
                periodic_incidence["daily"][TOTAL_SEVERE]   += 1;
            }

            if (infec->getCriticalTime() == now) {
                periodic_incidence["daily"][INTRO_CRITICAL] += intro;
                periodic_incidence["daily"][TOTAL_CRITICAL] += 1;
            }

            if (infec->getDeathTime() == now) {
                periodic_incidence["daily"][INTRO_DEATH]    += intro;
                periodic_incidence["daily"][TOTAL_DEATH]    += 1;
            }

            // Prevalence
            periodic_prevalence[INTRO_INF_PREV]      += intro and infec->isInfected(now);
            periodic_prevalence[INTRO_CASE_PREV]     += intro and infec->isSymptomatic(now);
            periodic_prevalence[INTRO_SEVERE_PREV]   += intro and infec->isSevere(now);
            periodic_prevalence[INTRO_CRITICAL_PREV] += intro and infec->isCritical(now);
            periodic_prevalence[INTRO_DEATH_PREV]    += intro and infec->isDead(now);

            periodic_prevalence[TOTAL_INF_PREV]      += infec->isInfected(now);
            periodic_prevalence[TOTAL_CASE_PREV]     += infec->isSymptomatic(now);
            periodic_prevalence[TOTAL_SEVERE_PREV]   += infec->isSevere(now);
            periodic_prevalence[TOTAL_CRITICAL_PREV] += infec->isCritical(now);
            periodic_prevalence[TOTAL_DEATH_PREV]    += infec->isDead(now);
        }
    }

    periodic_output(par, periodic_incidence, periodic_prevalence, date, process_id, epi_sizes);
    return;
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
Date,rcase,rdeath,rhosp
2020-03-02,2,0,
2021-04-29,5394,81,
*/
map<size_t, int> parse_emp_data_file(const Parameters* par, const string emp_data_file) {
    vector< vector<string> > emp_data = read_2D_vector_file(emp_data_file, ',');
    map<size_t, int> recast_emp_data;

    for (vector<string> &v : emp_data) {
        if (v[0] == "Date") { continue; }
        recast_emp_data[Date::to_sim_day(par->startJulianYear, par->startDayOfYear, v[0])] = stoi(v[1]);
    }
    return recast_emp_data;
}

enum FittingDataSource { EMP, SIM, NUM_OF_FITTING_DATA_SOURCES };
pair<bool, vector<double> > check_social_contact_fit(map<size_t, int> &emp_data, const vector<size_t> &sim_rcases, const size_t first_sim_day_of_month, const size_t last_sim_day_of_month,
                                                        const double fitting_threshold) {
    const size_t emp_start = emp_data.begin()->first;
    const size_t emp_end   = emp_data.end()->first;
    // if the emp data does not contain vlaues for the sim window, continue the simulation
    if ( ((first_sim_day_of_month < emp_start) and (last_sim_day_of_month < emp_start)) or
         ((first_sim_day_of_month > emp_end) and (last_sim_day_of_month > emp_end)) ) { return pair<bool, vector<double> >{ true, {-1.0, -1.0} }; }
    else { // check fit here: return true if fit is acceptable and false to re-simulate
        const double fl_per_10k = 1e4/20609673;//21538187; // TODO(alex): i dont like having these numbers here, but dont know how to solve it yet; will come back to this later
        const double sim_per_10k = 1e4/375474;

        const size_t overlap_start = max(emp_start, first_sim_day_of_month);
        const size_t overlap_end   = min(emp_end, last_sim_day_of_month);

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
            for (size_t day = 0; day < relevant_data[i].size(); ++day) {
                cumsum += relevant_data[i].at(day);
                cumsum_data[i].at(day) = cumsum;
            }
        }

//        cerr << "DAILY RCASES" << endl;
//        for (const vector<double> &v : relevant_data) {
//            for(const double &val : v) { cerr << val << ' '; }
//            cerr << endl;
//        }
//        cerr << endl;

//        cerr << "DAILY CUMUL RCASES" << endl;
//        cerr << "WIN DAY" << "\t\t" << "EMP" << "\t\t" << "SIM" << "\t\t" << "SIM-EMP" << "\t\t" << "SQR DIST" << endl;
        vector<double> sqr_dist;
        vector<double> dist;
        for (size_t i = overlap_start; i <= overlap_end; ++i) {
//            cerr << i << "\t\t" << cumsum_data[EMP].at(i) << "\t\t" << cumsum_data[SIM].at(i) << "\t\t";
            const double diff = cumsum_data[SIM].at(i) - cumsum_data[EMP].at(i);
//            cerr << diff << "\t\t" << (diff * diff) << endl;
            dist.push_back(diff);
            sqr_dist.push_back( (diff * diff) );
        }
//        cerr << endl;

        const double avg_sqr_dist = mean(sqr_dist);
        const double avg_dist = mean(dist);
        cerr << "AVG SQR DIST: " << avg_sqr_dist << " AVG DIST: " << avg_dist << endl << endl;
        //return pair<bool, vector<double> >{ true, {avg_dist, avg_sqr_dist} }; // use for testing/debugging
        if (avg_sqr_dist < fitting_threshold) {
            cerr << "ACCEPTABLE FIT" << endl;
            return pair<bool, vector<double> >{ true, {avg_dist, avg_sqr_dist} };
        } else {
            cerr << "UNACCEPTABLE FIT" << endl;
            return pair<bool, vector<double> >{ false, {avg_dist, avg_sqr_dist} };
        }

//        cerr << "TRAILING AVGS" << endl;
//        vector< vector<double> > trail_avgs(relevant_data.size(), vector<double>());
//        for (size_t i = 0; i < relevant_data.size(); ++i) { trail_avgs[i] = calc_trailing_avg(relevant_data[i], 7); }
//
//        for (const vector<double> &v : trail_avgs) {
//            for (const double &val : v) { cerr << val << ' '; }
//            cerr << endl;
//        }
//        cerr << endl;
//
//        cerr << "CENTERED AVGS" << endl;
//        vector< vector<double> > cen_avgs(relevant_data.size(), vector<double>());
//        for (size_t i = 0; i < relevant_data.size(); ++i) { cen_avgs[i] = calc_centered_avg(relevant_data[i], 7); }
//
//        for (const vector<double> &v : cen_avgs) {
//            for (const double &val : v) { cerr << val << ' '; }
//            cerr << endl;
//        }
//        cerr << endl;

//        cerr << "WEEKLY SUMS" << endl;
//        for (const vector<double> &v : relevant_data) {
//            double wk_sum = 0.0;
//            for (size_t i = 0; i < v.size(); ++i){
//                if ((i+1) % 7 == 0) {
//                    wk_sum += v[i];
//                    cerr << wk_sum << ' ';
//                    wk_sum = 0;
//                } else {
//                    wk_sum += v[i];
//                }
//            }
//            cerr << wk_sum << endl;
//        }
    }
}

vector<string> simulate_epidemic(const Parameters* par, Community* &community, const string process_id, const vector<string> mutant_intro_dates) {
    vector<int> epi_sizes;
    vector<string> daily_output_buffer;

    map<string, vector<int> > periodic_incidence = construct_tally();
    vector<int> periodic_prevalence(NUM_OF_PREVALENCE_REPORTING_TYPES, 0);
    vector<double> trailing_averages(par->runLength);
    const double pop_at_risk = min(community->getNumPeople(), par->numSurveilledPeople);

    vector<string> plot_log_buffer = {"date,sd,seasonality,vocprev1,vocprev2,cinf,closed,rcase,rdeath,inf,rhosp,Rt"};

    Date* date = community->get_date();
    vector<StrainType> strains = {WILDTYPE};

    vector<int> epi_sizes_ckpt;
    vector<string> daily_output_buffer_ckpt;
    map<string, vector<int> > periodic_incidence_ckpt;
    vector<int> periodic_prevalence_ckpt;
    vector<string> plot_log_buffer_ckpt;
    vector<StrainType> strains_ckpt;
    Community* community_ckpt   = nullptr;
    gsl_rng* RNG_ckpt           = nullptr;
    gsl_rng* REPORTING_RNG_ckpt = nullptr;

    vector<TimeSeriesAnchorPoint> social_contact_anchors = { // TODO(alex): do something better than just copying this here (but for now its fine)
        {"2020-01-01", 0.0},
        {"2020-03-10", 0.20},
        {"2020-03-15", 0.8},
        {"2020-04-01", 0.7},
        {"2020-05-01", 0.6},
        {"2020-06-01", 0.05},
        {"2020-07-01", 0.05},
        {"2020-08-01", 0.6},
        {"2020-09-01", 0.3},
        {"2020-10-01", 0.15},
        {"2020-11-01", 0.0},
        {"2020-12-01", 0.0},
        {"2021-01-01", 0.1},
        {"2021-02-01", 0.3},
        {"2021-02-28", 0.1},
        {"2021-03-01", 0.0},
        {"2021-04-01", 0.0},
        {"2021-05-01", 0.0}
    };

//int checkpointCount = 0;
    map<size_t, int> emp_data = parse_emp_data_file(par, "rcasedeath-florida.csv");

    for (; date->day() < (signed) par->runLength; date->increment()) {
        const size_t sim_day = date->day();
        //update_vaccinations(par, community, date);
        community->tick();

        if ( mutant_intro_dates.size() ) {
            if (*date >= mutant_intro_dates[0] and *date < mutant_intro_dates[1]) {
                //strains = {WILDTYPE, B_1_1_7};
                strains = {B_1_1_7};
            } else if (*date >= mutant_intro_dates[1]) {
                strains = {B_1_617_2};
            }
        }

        seed_epidemic(par, community, date, strains);
        const vector<size_t> infections         = community->getNumNewlyInfected();
        const vector<size_t> all_reported_cases = community->getNumDetectedCasesReport();
        const size_t reported_cases             = all_reported_cases[sim_day];
        const double trailing_avg               = calc_trailing_avg(all_reported_cases, sim_day, 7); // <= 7-day trailing average
        const vector<size_t> rhosp              = community->getNumDetectedHospitalizations();
        const vector<size_t> severe_prev        = community->getNumSeverePrev();
        const double cinf                       = accumulate(infections.begin(), infections.begin()+sim_day+1, 0.0);
        const double cAR                        = cinf/pop_at_risk; // cumulative attack rate (I hate this term)
        const vector<size_t> rdeaths            = community->getNumDetectedDeaths();


        const size_t rc_ct = accumulate(all_reported_cases.begin(), all_reported_cases.begin()+sim_day+1, 0);
        if (date->dayOfMonth()==1) cerr << "        rep sday        date  infinc  cAR     rcases  rcta7  crcases  rdeath  crdeath  sevprev   crhosp  closed  socdist\n";
        cerr << right
             << setw(11) << process_id
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
             << "  "     << setprecision(2) << par->timedInterventions.at(SOCIAL_DISTANCING).at(sim_day)
             << endl;

        //date,sd,seasonality,vocprev,cinf,closed,rcase,rdeath,Rt
        stringstream ss;
        ss << date->to_string({"yyyy", "mm", "dd"}, "-") << ","
           << par->timedInterventions.at(SOCIAL_DISTANCING).at(sim_day) << ","
           << par->seasonality.at(date->julianDay()-1) << ","
           << (float) community->getNumNewInfections(B_1_1_7)[sim_day]/infections[sim_day] << ","
           << (float) community->getNumNewInfections(B_1_617_2)[sim_day]/infections[sim_day] << ","
           << cAR << ","
           << community->getTimedIntervention(NONESSENTIAL_BUSINESS_CLOSURE, sim_day)<< ","
           << reported_cases*1e4/pop_at_risk << ","
           << rdeaths[sim_day]*1e4/pop_at_risk << ","
           << infections[sim_day]*1e4/pop_at_risk << ","
           << rhosp[sim_day]*1e4/pop_at_risk;
        plot_log_buffer.push_back(ss.str());


// if checkpoint should be stored now, cache all the stuff
// else if checkpoint should be restored now, switch to cached version, delete the current version
/*  Things to do in order to implement checkpoint caching with automated fitting of social contact curve:
        - check if the sim day is a checkpoint (i.e. last day of the month)
        - if at a checkpointing day --> cache simulation state
        - continue simulation for up to next checkpoint (i.e. next last day of the month)
        - check fit of social contact curve (HOW???)
            - if fit is good --> replace cache --> continue
            - if fit is bad  --> alter this month's social contact parameter --> reload cache --> continue

        - pass ap into simulate_epidemic by reference
*/

        if (sim_day == 0) { // create initial cache
cerr << "INITIAL SIM CACHE" << endl;
            epi_sizes_ckpt           = epi_sizes;
            daily_output_buffer_ckpt = daily_output_buffer;
            periodic_incidence_ckpt  = periodic_incidence;
            periodic_prevalence_ckpt = periodic_prevalence;
            plot_log_buffer_ckpt     = plot_log_buffer;
            strains_ckpt             = strains;
            if (community_ckpt) { delete community_ckpt; }
            community_ckpt           = new Community(*community);
            if (RNG_ckpt)           { gsl_rng_free(RNG_ckpt); }
            if (REPORTING_RNG_ckpt) { gsl_rng_free(REPORTING_RNG_ckpt); }
            RNG_ckpt                 = gsl_rng_clone(RNG);
            REPORTING_RNG_ckpt       = gsl_rng_clone(REPORTING_RNG);
        } else if (date->endOfMonth()) {
            // if we are here, sim_day is the last day of the month (day of month - 1 get us to the start of the month)
            const size_t first_sim_day_of_month = ((signed long) sim_day) - ((signed long) date->dayOfMonth() - 1) > 0
                                                  ? ((signed long) sim_day) - ((signed long) date->dayOfMonth() - 1)
                                                  : 0;

            const double fitting_threshold = 1e2;
            pair<bool, vector<double> > social_contact_fit = check_social_contact_fit(emp_data, all_reported_cases, first_sim_day_of_month, sim_day, fitting_threshold);

            // calculate the change we want to make to the social contact parameter
            const double avg_dist     = get<vector<double>>(social_contact_fit)[0];
            const double avg_sqr_dist = get<vector<double>>(social_contact_fit)[1];

            // if avg_dist > 0 --> sim was higher than emp --> need to increase social contact param for that month (decreases transmission)
            // if avg_dist < 0 --> sim was lower than emp --> need to decrease social contact param for that month (increases transmission)
            // adjust parameter using logistic curve: y = 0.1/(1 + e^(-0.1(x - 50))) where 0.1 is the max we would change the parameter at one time and 50 is half of the fitting threshold (100 rn)
            const double max_adj = 0.1;
            double social_contact_param_adj = (max_adj)/(1 + exp(-max_adj * (avg_sqr_dist - (fitting_threshold/2))));
            social_contact_param_adj = (avg_dist > 0) ? social_contact_param_adj : -social_contact_param_adj;

            // TODO(alex): auto fitting needs to be able to adjust transmissibility when fit is bad and social contact param is already 0
            if (get<bool>(social_contact_fit) or
                ((get<bool>(social_contact_fit) == false) and (community->social_distancing(first_sim_day_of_month) == 0)) ) { // replace cache
cerr << "REPLACING CACHED SIM" << endl;
                epi_sizes_ckpt           = epi_sizes;
                daily_output_buffer_ckpt = daily_output_buffer;
                periodic_incidence_ckpt  = periodic_incidence;
                periodic_prevalence_ckpt = periodic_prevalence;
                plot_log_buffer_ckpt     = plot_log_buffer;
                strains_ckpt             = strains;
                if (community_ckpt) { delete community_ckpt; }
                community_ckpt           = new Community(*community);
                if (RNG_ckpt)           { gsl_rng_free(RNG_ckpt); }
                if (REPORTING_RNG_ckpt) { gsl_rng_free(REPORTING_RNG_ckpt); }
                RNG_ckpt                 = gsl_rng_clone(RNG);
                REPORTING_RNG_ckpt       = gsl_rng_clone(REPORTING_RNG);
            } else { // change this month's social contact parameter and restore cache
                // iterate through ap to find points with date where: first day of current month <= date <= last day of current month
                for (TimeSeriesAnchorPoint &tsap : social_contact_anchors) {
                    const size_t tsap_sim_day = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, tsap.date);
                    if ( (tsap_sim_day >= first_sim_day_of_month) and (tsap_sim_day <= sim_day) ) { tsap.value += social_contact_param_adj; }
                }
                // restore cache
cerr << "RESTORING" << endl;
                epi_sizes                = epi_sizes_ckpt;
                daily_output_buffer      = daily_output_buffer_ckpt;
                periodic_incidence       = periodic_incidence_ckpt;
                periodic_prevalence      = periodic_prevalence_ckpt;
                plot_log_buffer          = plot_log_buffer_ckpt;
                strains                  = strains_ckpt;
                if (community) { delete community; }
                community                = new Community(*community_ckpt);
                date                     = community->get_date();
                if (RNG_ckpt) { gsl_rng_memcpy(RNG, RNG_ckpt); }
                if (REPORTING_RNG_ckpt) { gsl_rng_memcpy(REPORTING_RNG, REPORTING_RNG_ckpt); }

                // re-define social contact curve
                community->setSocialDistancingTimedIntervention(social_contact_anchors);
            }
        }

/*        if (sim_day == 30) { // checkpoint now
cerr << "CHECKPOINTING" << endl;
            epi_sizes_ckpt           = epi_sizes;
            daily_output_buffer_ckpt = daily_output_buffer;
            periodic_incidence_ckpt  = periodic_incidence;
            periodic_prevalence_ckpt = periodic_prevalence;
            plot_log_buffer_ckpt     = plot_log_buffer;
            strains_ckpt             = strains;
            if (community_ckpt) { delete community_ckpt; }
            community_ckpt           = new Community(*community);
            RNG_ckpt                 = gsl_rng_clone(RNG);
            REPORTING_RNG_ckpt       = gsl_rng_clone(REPORTING_RNG);
        } else if (sim_day == 45 and checkpointCount++ < 1) { // reset now
cerr << "RESTORING" << endl;
            epi_sizes                = epi_sizes_ckpt;
            daily_output_buffer      = daily_output_buffer_ckpt;
            periodic_incidence       = periodic_incidence_ckpt;
            periodic_prevalence      = periodic_prevalence_ckpt;
            plot_log_buffer          = plot_log_buffer_ckpt;
            strains                  = strains_ckpt;
            if (community) { delete community; }
            community                = new Community(*community_ckpt);
            date                     = community->get_date();
            if (RNG_ckpt) { gsl_rng_memcpy(RNG, RNG_ckpt); }
            if (REPORTING_RNG_ckpt) { gsl_rng_memcpy(REPORTING_RNG, REPORTING_RNG_ckpt); }
        }*/

    }

    if (community_ckpt)     { delete community_ckpt; }
    if (RNG_ckpt)           { gsl_rng_free(RNG_ckpt); }
    if (REPORTING_RNG_ckpt) { gsl_rng_free(REPORTING_RNG_ckpt); }
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
    cerr << "icu deaths, total deaths, ratio: " << cdeath_icu << ", " << cdeath2 << ", " << cdeath_icu/cdeath2 << endl;

    cerr << "FINAL SOCIAL CONTACT PARAMS" << endl;
    cerr << "DATE\t\tVALUE" << endl;
    for (const TimeSeriesAnchorPoint &tsap : social_contact_anchors) { cerr << tsap.date << "\t\t" << tsap.value << endl; }

//  write_daily_buffer(plot_log_buffer, process_id, "plot_log.csv");
    //return epi_sizes;
    if (RNG) { gsl_rng_free(RNG); }
    if (REPORTING_RNG) { gsl_rng_free(REPORTING_RNG); }
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
