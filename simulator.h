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


const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);

// Predeclare local functions
Community* build_community(const Parameters* par);
void seed_epidemic(const Parameters* par, Community* community);
vector<int> simulate_epidemic(const Parameters* par, Community* community, const string process_id = "0");
void write_immunity_file(const Community* community, const string label, string filename, int runLength);
void write_immunity_by_age_file(const Community* community, const int year, string filename="");
void write_output(const Parameters* par, Community* community, vector<int> initial_susceptibles);
void write_daily_buffer( vector<string>& buffer, const string process_id, string filename);

Community* build_community(const Parameters* par) {
    Community* community = new Community(par);
    Person::setPar(par);
cerr << "Reading locations ... ";
    if (!community->loadLocations(par->locationFilename)) {
        cerr << "ERROR: Could not load locations" << endl;
        exit(-1);
    }
cerr << "done.\n";
cerr << "Reading population ... ";
    if (!community->loadPopulation(par->populationFilename)) {
        cerr << "ERROR: Could not load population" << endl;
        exit(-1);
    }

//cerr << "done.\n  Now sleeping for 20s so ram usage can be checked.\n";
//sleep(20);
    if (!par->abcVerbose) {
        cerr << community->getNumPeople() << " people" << endl;
    }

    /*if (!par->bSecondaryTransmission) {
        community->setNoSecondaryTransmission();
    }*/

    return community;
}


void seed_epidemic(const Parameters* par, Community* community) {
    // epidemic may be seeded with initial exposure OR initial infection
    bool attempt_initial_infection = true;
    // Normal usage, to simulate epidemic
    if (par->numInitialExposed > 0) {
        attempt_initial_infection = false;
        for (size_t i=0; i<par->numInitialExposed; i++)
            community->infect(gsl_rng_uniform_int(RNG, community->getNumPeople()), 0);
    }
    if (attempt_initial_infection) {
        // Useful for estimating R0
        if(par->numInitialInfected > 0) {
            int count = community->getNumInfected(0);

            // must infect initialInfected persons -- this bit is mysterious
            while (community->getNumInfected(0) < count + par->numInitialInfected) {
                community->infect(gsl_rng_uniform_int(RNG, community->getNumPeople()), 0);
            }
        }
    }
    return;
}

/*
void write_yearly_people_file(const Parameters* par, const Community* community, int time) {
    ofstream yearlyPeopleOutputFile;
    ostringstream ssFilename;
    ssFilename << par->yearlyPeopleOutputFilename << ((int)(time/365)) << ".csv";
    cerr << "outputing yearly people information to " << ssFilename.str() << endl;
    yearlyPeopleOutputFile.open(ssFilename.str().c_str());
    if(yearlyPeopleOutputFile.fail()) {
        cerr << "ERROR: People file '" << par->yearlyPeopleOutputFilename << "' cannot be open for writing." << endl;
        exit(-1);
    }
    yearlyPeopleOutputFile << "pid,infectiontime,symptomtime,withdrawtime,recoverytime,immdenv1,immdenv2,immdenv3,immdenv4" << endl;
    for (Person* p: community->getPeople()) {
        for (int j=p->getNumNaturalInfections()-1; j>=0; j--) {
            yearlyPeopleOutputFile << p->getID() << ","
                << p->getInfectedTime(j) << ","
                << p->getSymptomTime(j) << ","
//                << p->getWithdrawnTime(j) << ","
//                << p->getRecoveryTime(j) << ",";
                yearlyPeopleOutputFile << (p->isNaive()?0:1) << endl;
        }
    }
    yearlyPeopleOutputFile.close();
    return;
}*/


void _aggregator(map<string, vector<int> >& periodic_incidence, string key) {
    for (unsigned int i = 0; i < periodic_incidence["daily"].size(); ++i) periodic_incidence[key][i] += periodic_incidence["daily"][i];
}


void _reporter(stringstream& ss, map<string, vector<int> > &periodic_incidence, vector<int> &periodic_prevalence, const Parameters* par, const string process_id, const string label, const int value, string key) {
        ss << process_id << dec << " " << par->serial << label << value << " ";
        for (auto v: periodic_incidence[key]) ss << v << " ";
        if(key=="daily") for (auto v: periodic_prevalence) ss << v << " ";
        for (auto v: par->reportedFraction) ss << v << " ";
}


void periodic_output(const Parameters* par, map<string, vector<int> > &periodic_incidence, vector<int> &periodic_prevalence, const Date& date, const string process_id, vector<int>& epi_sizes) {
    stringstream ss;
//if (date.day() >= 25*365 and date.day() < 36*365) {
//if (date.day() >= 116*365) {                         // daily output starting in 1995, assuming Jan 1, 1879 simulation start
//if (date.day() >= 99*365 and date.day() < 105*365) { // daily output for summer/winter IRS comparison
    if (par->dailyOutput) {
        _reporter(ss, periodic_incidence, periodic_prevalence, par, process_id, " day: ", date.day(), "daily"); ss << endl;
    }
//}
    vector<int> dummy;
     if (par->periodicOutput) {
        _aggregator(periodic_incidence, "n_day");
        const int n = par->periodicOutputInterval;
        if (date.endOfPeriod(n)) {
            _reporter(ss, periodic_incidence, dummy, par, process_id, " " + to_string(n) + "_day: ", date.nDayPeriod(n), "n_day"); ss << endl;
            periodic_incidence["n_day"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
        }
    }

    if (par->weeklyOutput) {
        _aggregator(periodic_incidence, "weekly");
        if (date.endOfWeek()) {
            _reporter(ss, periodic_incidence, dummy, par, process_id, " week: ", date.week(), "weekly"); ss << endl;
            periodic_incidence["weekly"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
        }
    }

    if (par->monthlyOutput) {
        _aggregator(periodic_incidence, "monthly");
        if (date.endOfMonth()) {
            _reporter(ss, periodic_incidence, dummy, par, process_id, " month: ", date.julianMonth(), "monthly"); ss << endl;
            periodic_incidence["monthly"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
        }
    }

    // handle several things that happen yearly
    _aggregator(periodic_incidence, "yearly");
    if (date.endOfYear()) {
        if (par->abcVerbose) {
            cout << process_id << dec << " " << par->serial << " T: " << date.day() << " annual: ";
            for (auto v: periodic_incidence["yearly"]) { cout << v << " "; } cout << endl;
        }

        epi_sizes.push_back(periodic_incidence["yearly"][2]);

        if (par->yearlyOutput) { _reporter(ss, periodic_incidence, dummy, par, process_id, " year: ", date.year(), "yearly"); ss << endl; }
        periodic_incidence["yearly"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
    }

    periodic_incidence["daily"] = vector<int>(NUM_OF_INCIDENCE_REPORTING_TYPES, 0);
    string output = ss.str();
    //fputs(output.c_str(), stderr);
    fputs(output.c_str(), stdout);
}

void update_vaccinations(const Parameters* par, Community* community, const Date &date) {
    const int doseInterval = par->vaccineDoseInterval;
    assert(doseInterval > 0); // neg is nonsensical, 0 is disallowed due to mod operation
    //const int boostInterval = par->vaccineBoostingInterval;
    for (CatchupVaccinationEvent cve: par->catchupVaccinationEvents) {
        // Normal, initial vaccination -- boosting, multiple doses handled in Community::tick()
        if (date.day() == cve.simDay) {
            if (not par->abcVerbose) cerr << "vaccinating " << cve.coverage*100 << "% of age " << cve.age << " on day " << cve.simDay << endl;
            community->vaccinate(cve);
        }
    }
}


int seed_epidemic(const Parameters* par, Community* community, const Date &date) {
    int introduced_infection_ct = 0;
    const int numperson = community->getNumPeople();
    const int ai_year_lookup = date.year() % par->annualIntroductions.size();
    const double intros = par->annualIntroductions[ai_year_lookup];
    //const int de_year_lookup = date.year() % par->numDailyExposed.size();
    //const double weight = par->dailyExposed[de_year_lookup];
    const double annual_intros_weight = par->annualIntroductionsCoef;
    const double expected_num_exposed = /*weight * */ annual_intros_weight * intros;
    if (expected_num_exposed > 0) {
        assert(expected_num_exposed <= numperson);
        const int num_exposed = gsl_ran_poisson(RNG, expected_num_exposed);
        for (int i=0; i<num_exposed; i++) {
            // gsl_rng_uniform_int returns on [0, numperson-1]
            int transmit_to_id = gsl_rng_uniform_int(RNG, numperson);
            if (community->infect(transmit_to_id, date.day())) {
                introduced_infection_ct++;
            }
        }
    }
    return introduced_infection_ct;
}


void advance_simulator(const Parameters* par, Community* community, Date &date, const string process_id, map<string, vector<int> > &periodic_incidence, vector<int> &periodic_prevalence, vector<int> &epi_sizes) {
    community->tick(date.day());

    seed_epidemic(par, community, date);

    for (Person* p: community->getPeople()) {
        const int now = date.day();
        if (p->isInfected(now) or p->isSymptomatic(now)) {
            const Infection* infec = p->getInfection();
            bool intro  = not infec->isLocallyAcquired();
            bool symp   = infec->isSymptomatic(now);      // subset of all infections
            bool severe = infec->isSevere(now);           // subset of symptomatic
            bool crit   = infec->isCritical(now);         // subset of severe
            bool dead   = infec->isDead(now);             // exclusive of other states

            // Prevalence
            periodic_prevalence[INTRO_INF_PREV]      += intro;
            periodic_prevalence[INTRO_CASE_PREV]     += intro and symp;
            periodic_prevalence[INTRO_SEVERE_PREV]   += intro and severe;
            periodic_prevalence[INTRO_CRITICAL_PREV] += intro and crit;
            periodic_prevalence[INTRO_DEATH_PREV]    += intro and dead;

            periodic_prevalence[TOTAL_INF_PREV]      += 1; // if you're here, then you're infected
            periodic_prevalence[TOTAL_CASE_PREV]     += symp;
            periodic_prevalence[TOTAL_SEVERE_PREV]   += severe;
            periodic_prevalence[TOTAL_CRITICAL_PREV] += crit;
            periodic_prevalence[TOTAL_DEATH_PREV]    += dead;

            // Incidence
            if (p->isNewlyInfected(date.day())) {
                periodic_incidence["daily"][INTRO_INF]      += intro;
                periodic_incidence["daily"][INTRO_CASE]     += intro and symp;
                periodic_incidence["daily"][INTRO_SEVERE]   += intro and severe;
                periodic_incidence["daily"][INTRO_CRITICAL] += intro and crit;
                periodic_incidence["daily"][INTRO_DEATH]    += intro and dead;

                periodic_incidence["daily"][TOTAL_INF]      += 1;
                periodic_incidence["daily"][TOTAL_CASE]     += symp;
                periodic_incidence["daily"][TOTAL_SEVERE]   += severe;
                periodic_incidence["daily"][TOTAL_CRITICAL] += crit;
                periodic_incidence["daily"][TOTAL_DEATH]    += dead;
            }
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


vector<int> simulate_epidemic(const Parameters* par, Community* community, const string process_id) {
    vector<int> epi_sizes;
    Date date(par);

    vector<string> daily_output_buffer;

/*    if (not par->abcVerbose) {
        daily_output_buffer.push_back("day,year,id,age,location,vaccinated,symptomatic,severity");
    }*/

    map<string, vector<int> > periodic_incidence = construct_tally();
    vector<int> periodic_prevalence(NUM_OF_PREVALENCE_REPORTING_TYPES, 0);

    for (; date.day() < par->runLength; date.increment()) {
        update_vaccinations(par, community, date);
        advance_simulator(par, community, date, process_id, periodic_incidence, periodic_prevalence, epi_sizes);
/*        if ((date.julianDay() == ((sero_prev_aggregation_julian_start+364) % 365 ) + 1)) { // +1 because julianDay is [1,365])), avg(avg(interventions are specified on [0,364]
            // tally current seroprevalence stats
            int vaccinated_tally = 0;
            vector<double> inf_ct_tally(5,0.0);         // [0,4] past infections
            for (Person* p: community->getPeople()) {
                // ignoring possible seroconversion due to vaccine
                vaccinated_tally += p->isVaccinated();
                const int ct = p->getNumNaturalInfections();
                ++inf_ct_tally[ct];
            }
            double N = community->getNumPeople();
            cerr << "vaccinated N, fraction: " << vaccinated_tally << " " << (double) vaccinated_tally / N << endl;
        }*/
    }
    return epi_sizes;
}


bool fileExists(const std::string& filename) {
    struct stat buf;
    return stat(filename.c_str(), &buf) != -1;
}


void write_daily_buffer( vector<string>& buffer, const string process_id, string filename = "" ) {
    if (filename == "") {
        stringstream ss_filename;
        ss_filename << "daily_output." << process_id;
        filename = ss_filename.str();
    }

    string all_output;
    for (const auto &line : buffer) all_output += (line + "\n");

    if (fileExists(filename)) {
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


void write_output(const Parameters* /*par*/, Community* /*community*/, vector<int> /*numInitialSusceptible*/) {
/*    if (!par->bSecondaryTransmission) {
        // outputs
        //   number of households infected
        //   age of index case
        //   age(s) of secondary cases
        vector<int> numCurrentSusceptible = community->getNumSusceptible();
        cout << (numInitialSusceptible - numCurrentSusceptible) << " ";
        //    cout << "secondary infections" << endl;

        int ages[100];
        int times[100];
        int numages=0;
        int indexage=-1;
        int homeids[100];
        int numhomes=0;
        // ages of infected people
        for (Person* p: community->getPeople()) {
            int t = p->getInfectedTime();
            if (t>=0) {
                if (t==0)
                    indexage = p->getAge();
                else {
                    ages[numages] = p->getAge();
                    times[numages] = t;
                    numages++;
                }
                int homeid = p->getHomeID();
                bool bFound = false;
                for (int j=0; j<numhomes; j++)
                    if (homeids[j]==homeid)
                        bFound = true;
                if (!bFound)
                    homeids[numhomes++] = homeid;
            }
        }
        cout << indexage << " " << numhomes << " " << numages;
        for (int i=0; i<numages; i++)
            cout << " " << ages[i];
        for (int i=0; i<numages; i++)
            cout << " " << times[i];
        cout << endl;
    }

    // output daily infected/symptomatic file
    if (par->dailyOutputFilename.length()>0) {
        cerr << "outputing daily infected/symptomatic information to " << par->dailyOutputFilename << endl;
        ofstream dailyOutputFile;
        dailyOutputFile.open(par->dailyOutputFilename.c_str());
        if(dailyOutputFile.fail()) {
            cerr << "ERROR: Daily file '" << par->dailyOutputFilename << "' cannot be open for writing." << endl;
            exit(-1);
        }
        dailyOutputFile << "day,newly infected DENV1,newly infected DENV2,newly infected DENV3,newly infected DENV4,"
                  << "newly symptomatic DENV1,newly symptomatic DENV2,newly symptomatic DENV3,newly symptomatic DENV4" << endl;
        vector< vector<int> > infected =    community->getNumNewlyInfected();
        vector< vector<int> > symptomatic = community->getNumNewlySymptomatic();
        for (int t=0; t<par->runLength; t++) {
            dailyOutputFile << t << ",";
            for (int i=0; i<NUM_OF_SEROTYPES; i++)   dailyOutputFile << infected[i][t] << ",";
            for (int i=0; i<NUM_OF_SEROTYPES-1; i++) dailyOutputFile << symptomatic[i][t] << ",";
            dailyOutputFile << symptomatic[NUM_OF_SEROTYPES-1][t] << endl;
        }
        dailyOutputFile.close();
    }

    // output people file
    if (par->peopleOutputFilename.length()>0) {
        cerr << "outputing people information to " << par->peopleOutputFilename << endl;
        ofstream peopleOutputFile;
        peopleOutputFile.open(par->peopleOutputFilename.c_str());
        if(peopleOutputFile.fail()) {
            cerr << "ERROR: People file '" << par->peopleOutputFilename << "' cannot be open for writing." << endl;
            exit(-1);
        }
        peopleOutputFile << "pid,serotype,infectiontime,symptomtime,withdrawtime,recoverytime,immdenv1,immdenv2,immdenv3,immdenv4,vaccinated" << endl;
        for (Person* p: community->getPeople()) {
            for (int j=p->getNumNaturalInfections()-1; j>=0; j--) {
                peopleOutputFile << p->getID() << ","
                    << 1 + (int) p->getSerotype(j) << ","
                    << p->getInfectedTime(j) << ","
                    << p->getSymptomTime(j) << ","
                    << p->getWithdrawnTime(j) << ","
                    << p->getRecoveryTime(j) << ",";
                for (int s = 0; s < NUM_OF_SEROTYPES; ++s) {
                    peopleOutputFile << (p->isSusceptible((Serotype) s)?0:1) << ",";
                }
                    peopleOutputFile << (p->isVaccinated()?1:0) << endl;
            }
        }
        peopleOutputFile.close();
    }*/
}
