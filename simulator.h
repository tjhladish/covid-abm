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
#include "behavior_tuner.h"

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
gsl_rng* VAX_RNG = gsl_rng_alloc(gsl_rng_mt19937);

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
        // cerr << "pop size, sampled size, infected size: " << pop_size << ", " << k << ", " << inf_ct << endl;
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
    ledger->plot_log_buffer = {"serial,date,sd,seasonality,vocprev1,vocprev2,vocprev3,cinf,closed,rcase,rdeath,inf,rhosp,VES,brkthruRatio,vaxInfs,unvaxInfs,hospInc,hospPrev,icuInc,icuPrev,vaxHosp,unvaxHosp,std_doses,urg_doses,cov1,cov2,cov3,seroprev,symp_infs,sevr_infs,crit_infs,all_deaths,Rt"};
    //ledger->plot_log_buffer = {"date,sd,seasonality,vocprev1,vocprev2,cinf,closed,rcase,rdeath,inf,rhosp,Rt"};

    Date* date = community->get_date();
    ledger->strains = {50.0, 0.0, 0.0, 0.0}; // initially all WILDTYPE
    assert(ledger->strains.size() == NUM_OF_STRAIN_TYPES);

    if (par->behavioral_autotuning) {
        // create tuner and initialize first simulation cache
        tuner = initialize_behavior_autotuning(par);
        sim_cache = new SimulationCache(community, ledger, RNG, REPORTING_RNG, VAX_RNG);    // TODO: change to quick cache
        first_tuning_window_setup(par, community, tuner, social_distancing_anchors);
    } else if (not par->behaviorInputFilename.empty()) {
        // filename provided for a dataset with behavioral values to use
        init_behavioral_vals_from_file(par, community);
    }

    bool restore_occurred = false; // relevant for behavior autotuning
    for (; date->day() < (signed) par->runLength; date->increment()) {
        if (par->behavioral_autotuning) {
           behavior_autotuning(par, community, date, ledger, tuner, sim_cache, social_distancing_anchors, restore_occurred);
        }
        const size_t sim_day = date->day();
if (sim_day == 0) { seed_epidemic(par, community, WILDTYPE); }
//if (*date == "2021-12-01") { gsl_rng_set(RNG, par->randomseed); // use something like this if we want to only have dynamic uncertainty beyond some date

        community->tick();


        if ( mutant_intro_dates.size() ) {
            if (*date >= mutant_intro_dates[0] and *date < mutant_intro_dates[1]) {
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
        const vector<size_t> symp_infs          = community->getNumNewlySymptomatic();
        const vector<size_t> sevr_infs          = community->getNumNewlySevere();
        const vector<size_t> crit_infs          = community->getNumNewlyCritical();
        const vector<size_t> deaths             = community->getNumNewlyDead();
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

        Vac_Campaign* vc    = community->getVac_Campaign() ? community->getVac_Campaign() : nullptr;
        const int std_doses = vc ? vc->get_std_doses_used(sim_day) : 0;
        const int urg_doses = vc ? vc->get_urg_doses_used(sim_day) : 0;
        const int all_doses = vc ? vc->get_all_doses_used(sim_day) : 0;

        const size_t rc_ct = accumulate(all_reported_cases.begin(), all_reported_cases.begin()+sim_day+1, 0);
        map<string, double> VE_data = community->calculate_vax_stats(sim_day);
        const double seroprev = community->doSerosurvey(NATURAL, community->getPeople(), sim_day);
        //vector<string> inf_by_loc_keys = {"home", "social", "work_staff", "patron", "school_staff", "student", "hcw", "patient", "ltcf_staff", "ltcf_resident"};
        //for (string key : inf_by_loc_keys) {
        //    cerr << "infLoc " << date->to_ymd() << ' ' << key << ' ' << community->getNumNewInfectionsByLoc(key)[sim_day] << endl;
        //}
        if (not par->behavioral_autotuning) {
            if (date->dayOfMonth()==1) cerr << "        rep sday        date  infinc  cAR     rcases  rcta7  crcases  rdeath  crdeath  sevprev   crhosp  closed  socdist  cov_1  cov_2  cov_3  std_doses  urg_doses  all_doses  quar%\n";
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
                << "  "     << setw(7) <<          setprecision(2) << left << community->getTimedIntervention(SOCIAL_DISTANCING, sim_day) << right
                << "  "     << setw(5) << fixed << setprecision(2) << left << VE_data["dose_1"] << right << defaultfloat
                << "  "     << setw(5) << fixed << setprecision(2) << left << VE_data["dose_2"] << right << defaultfloat
                << "  "     << setw(5) << fixed << setprecision(2) << left << VE_data["dose_3"] << right << defaultfloat
                << setw(11) << std_doses
                << setw(11) << urg_doses
                << setw(11) << all_doses
                << "  "     << setw(7) << fixed << setprecision(2) << left << community->getNumPeopleQuarantining(sim_day) << right << defaultfloat
                << setprecision(6)
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
           << VE_data["unvaxHosp"]*1e4/pop_at_risk << ","
           << std_doses*1e4/pop_at_risk << ","
           << urg_doses*1e4/pop_at_risk << ","
           << VE_data["dose_1"] << ","
           << VE_data["dose_2"] << ","
           << VE_data["dose_3"] << ","
           << seroprev << ","
           << symp_infs*1e4/pop_at_risk << ","
           << sevr_infs*1e4/pop_at_risk << ","
           << crit_infs*1e4/pop_at_risk << ","
           << deaths*1e4/pop_at_risk;
        ledger->plot_log_buffer.push_back(ss.str());
    }

    const vector<size_t> sim_reported_cases = community->getNumDetectedCasesReport();
    if (par->behavioral_autotuning) {
        write_anchors_to_file(par, social_distancing_anchors);
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

cerr_vector(community->getTimedIntervention(SOCIAL_DISTANCING)); cerr << endl;

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


void import_csv_to_db(string filename, string table, string db) {
    stringstream ss;
    ss << "sqlite3 " << db << " '.mode csv' '.import " << filename << ' ' << table << "'";
    string cmd_str = ss.str();
    int retval = system(cmd_str.c_str());
    if (retval == -1) { cerr << "System failed to import " << table << " data to db\n"; }
    return;
}

void generate_sim_data_db(const Parameters* par, const Community* community, const unsigned long int serial, vector<string> tables) {
    vector<stringstream> filenames(tables.size());
    for (size_t i = 0; i < tables.size(); ++i) {
        filenames[i] << "./" << tables[i] << "_" << serial << ".csv";
    }

    map<string, ofstream> ofiles;
    for (size_t i = 0; i < tables.size(); ++i) {
        ofiles[tables[i]] = ofstream(filenames[i].str(), std::ios::trunc);
    }

    bool all_files_open = true;
    for (const auto& [table, ofile] : ofiles) {
        all_files_open = all_files_open and (bool)ofile;
    }
    if (not all_files_open) { cerr << "FILES FAILED TO OPEN" << endl; exit(-1); }

    for (Person* p : community->getPeople()) {
        for (Infection* inf : p->getInfectionHistory()) {
            if (not inf) { continue; }
            int inf_place_id = inf->getInfectedPlace() ? inf->getInfectedPlace()->getID() : -1;
            int inf_by_id    = inf->getInfectedBy() ? inf->getInfectedBy()->getID() : -1;
            int inf_owner_id = inf->getInfectionOwner() ? inf->getInfectionOwner()->getID() : -1;

            ofiles["infection_history"] << inf << ','
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
                                        << inf->secondary_infection_tally() << "\n";

            for (Infection* sec_inf : inf->get_infections_caused()) {
                ofiles["secondary_infections"] << inf << ',' << sec_inf << "\n";
            }

            if (inf->getDetection()) {
                Detection* inf_det = inf->getDetection();
                ofiles["infection_detection"] << inf_det << ','
                                              << inf << ','
                                              << inf_det->detected_state << ','
                                              << inf_det->reported_time << "\n";
            }
        }

        for (int vax = 0; vax < (int) p->getVaccinationHistory().size(); ++vax) {
            ofiles["vaccination_history"] << p->getID() << ','
                                          << community->getVac_Campaign()->get_age_bin(p->getAge()) << ','
                                          << vax << ','
                                          << p->getVaccinationHistory()[vax] << ','
                                          << Date::to_ymd(p->getVaccinationHistory()[vax], par) << "\n";
        }
    }

    map<int, int> bin_pops = community->getVac_Campaign()->get_unique_age_bin_pops();
    Dose_Ptrs std_doses = community->getVac_Campaign()->get_std_doses_available();
    Dose_Ptrs urg_doses = community->getVac_Campaign()->get_urg_doses_available();

    Dose_Vals std_doses_used = community->getVac_Campaign()->get_std_doses_used();
    Dose_Vals urg_doses_used = community->getVac_Campaign()->get_urg_doses_used();

    for (int bin : community->getVac_Campaign()->get_unique_age_bins()) {
        ofiles["age_bins"] << bin << ','
                           << bin_pops[bin] << "\n";

        for (int day = 0; day < (int) par->runLength; ++day) {
            for (int dose = 0; dose < par->numVaccineDoses; ++dose) {
                ofiles["doses_available"] << day << ','
                                          << Date::to_ymd(day , par) << ','
                                          << dose << ','
                                          << bin << ','
                                          << *std_doses[day][dose][bin] << ','
                                          << *urg_doses[day][dose][bin] << "\n";

                ofiles["doses_used"] << day << ','
                                     << Date::to_ymd(day , par) << ','
                                     << dose << ','
                                     << bin << ','
                                     << std_doses_used[day][dose][bin] << ','
                                     << urg_doses_used[day][dose][bin] << "\n";
            }
        }
    }

    for (auto& [table, ofile] : ofiles) { ofile.close(); }

    stringstream ss, db;
    db << "sim_data_" << serial << ".sqlite";
    ss << "sqlite3 " << db.str() << " '.read gen_sim_db.sql'";

    string cmd_str = ss.str();
    int retval = system(cmd_str.c_str());
    if (retval == -1) { cerr << "System call to `sqlite3 " << db.str() << " '.read gen_sim_db.sql'` failed\n"; }

    for (size_t i = 0; i < tables.size(); ++i) {
        import_csv_to_db(filenames[i].str(), tables[i], db.str());
    }

    ss.str(string());
    ss << "rm";
    for (size_t i = 0; i < filenames.size(); ++i) { ss << ' ' << filenames[i].str(); }

    cmd_str = ss.str();
    retval = system(cmd_str.c_str());
    if (retval == -1) { cerr << "System call to delete infection csv files failed\n"; }
}
