#include <cstdlib>
#include <cstring>
#include <climits>
#include <iostream>
#include <string>
#include <iterator>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <math.h>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Person.h"
#include "Location.h"
#include "Community.h"
#include "Parameters.h"

using namespace covid::standard;
using covid::util::mean;
using covid::util::choice;

const Parameters* Community::_par;
Date* Community::_date;
vector< map<LocationType, map<Location*, int, Location::LocPtrComp>>> Community::_isHot;
set<Person*> Community::_revaccinate_set;
vector<size_t> Community::_numDetectedCasesOnset;
vector<size_t> Community::_numDetectedCasesReport;
vector<size_t> Community::_numDetectedDeaths;

//vector<Person*> Community::_peopleByAge;
//map<int, set<pair<Person*,Person*> > > Community::_delayedBirthdays;

int mod(int k, int n) { return ((k %= n) < 0) ? k+n : k; } // correct for non-negative n


Community::Community(const Parameters* parameters, Date* date) :
//    _exposedQueue(MAX_INCUBATION, vector<Person*>(0)),
    _numNewlyInfected(parameters->runLength), // +1 not needed; runLength is already a valid size
    _numNewlySymptomatic(parameters->runLength),
    _numNewlyDead(parameters->runLength),
    _numVaccinatedCases(parameters->runLength),
    _numSeverePrev(parameters->runLength),
    _numHospInc(parameters->runLength),
    _numHospPrev(parameters->runLength),
    _numIcuInc(parameters->runLength),
    _numIcuPrev(parameters->runLength)
    {
    _par = parameters;
    _date = date;
    _day = 0;
    //_mortality = NULL;
    //_bNoSecondaryTransmission = false;
    //_uniformSwap = true;
//    for (int a = 0; a<NUM_AGE_CLASSES; a++) _personAgeCohortSizes[a] = 0;
    _isHot.resize(_par->runLength);
    _numDetectedCasesOnset.resize(_par->runLength);
    _numDetectedCasesReport.resize(_par->runLength);
    _numDetectedDeaths.resize(_par->runLength);
    for (auto &e: _isHot) {
        for (size_t locType = 0; locType < NUM_OF_LOCATION_TYPES; ++locType) {
            e[(LocationType) locType] = {};
        }
    }
    timedInterventions = _par->timedInterventions;
}


Date* Community::get_date() {
    if (_date) {
        return _date;
    } else {
        cerr << "ERROR: Community::_date not defined in Community::get_date()\n";
        exit(-1);
    }
}


void Community::reset() { // used for r-zero calculations, to reset pop after a single intro
    // reset people
    for (Person* p: _people) {
/*        if (p->isWithdrawn(_day)) {
            p->getLocation(WORK_DAY)->addPerson(p,WORK_DAY);                    // goes back to work
            p->getLocation(HOME_MORNING)->removePerson(p,WORK_DAY);             // stops staying at home
        }*/
        p->resetImmunity(); // no past infections, not dead, not vaccinated
    }

    // reset locations
    for (auto &e: _isHot) e.clear();

    // clear community queues & tallies
//    _exposedQueue.clear();
    _numNewlyInfected.clear();
    _numNewlySymptomatic.clear();
    _numNewlyDead.clear();
    _numVaccinatedCases.clear();

//    _exposedQueue.resize(MAX_INCUBATION, vector<Person*>(0));
    _numNewlyInfected.resize(_par->runLength);
    _numNewlySymptomatic.resize(_par->runLength);
    _numNewlyDead.resize(_par->runLength);
    _numVaccinatedCases.resize(_par->runLength);
}


Community::~Community() {
    if (_date) delete _date;
    if (_people.size() > 0) { for (Person* p: _people) delete p; }

    Person::reset_ID_counter();
    _isHot.clear();

    for (unsigned int i = 0; i < _location.size(); i++ ) delete _location[i];
    _location.clear();
    Location::reset_ID_counter();

    for (auto& kv: _location_map) {
        kv.second.clear();
    }
    _location_map.clear();

//    _exposedQueue.clear();
    _personAgeCohort.clear();
    _numNewlyInfected.clear();
    _numNewlySymptomatic.clear();
    _numNewlyDead.clear();
    _numVaccinatedCases.clear();
}


bool Community::loadPopulation(string populationFilename, string comorbidityFilename, string immunityFilename) {
    ifstream iss(populationFilename);

    if (!iss) {
        cerr << "ERROR: " << populationFilename << " not found." << endl;
        return false;
    }
    string buffer;
//    int agecounts[NUM_AGE_CLASSES];
//    for (int i=0; i<NUM_AGE_CLASSES; i++) agecounts[i] = 0;

    istringstream line;
    // per IPUMS, expecting 1 for male, 2 for female for sex
    int id, hid, age, sex, did;//, empstat;
    while ( getline(iss,buffer) ) {
        line.clear();
        line.str(buffer);
        /*
        pid res_id sex age mov_id
        0 1948559 1 21 475292
        1 1948559 1 23 475292
        2 1948560 2 22 475292
        3 1948560 2 20 475292
        4 1948560 2 20 481010
        5 1948560 2 19 481010
        */
        if (line >> id >> hid >> sex >> age >> did ) { //>> empstat) {
            if (id != (signed) _people.size()) { // ensures indexing stays consistent
                cerr << "ERROR: Person ID's must be sequential integers starting at 0" << endl;
                return false;
            }
            Person* p = new Person();
            _people.push_back(p);
            p->setAge(age);
            p->setSex((SexType) sex);
            assert((signed) _location.size() > hid);
            if ((signed) _location.size() <= did) {
                cerr << line.str() << endl;
                exit(-1);
            }
            p->setHomeLoc(_location[hid]);
            _location[hid]->addPerson(p);
            if (did >= 0) {
                p->setDayLoc(_location[did]); // currently just any non-home daytime location--may be a school; -1 == NA
                _location[did]->addPerson(p);
            }
            if (_location[hid]->getType() == NURSINGHOME) p->setLongTermCare(true);
            //assert(age<NUM_AGE_CLASSES);
 //           agecounts[age]++;
        }
    }
    iss.close();

    //_peopleByAge = _people;
    //sort(_peopleByAge.begin(), _peopleByAge.end(), PerPtrComp());

    if (comorbidityFilename.length()>0) {
        iss.open(comorbidityFilename);

        if (!iss) {
            cerr << "ERROR: " << comorbidityFilename << " not found." << endl;
            return false;
        }

        bool com;
        while ( getline(iss,buffer) ) {
            line.clear();
            line.str(buffer);
            /*
            pid sex age undlycond
            0 1 33 0
            1 2 27 0
            2 2 58 0
            3 1 63 0
            4 2 49 1
            */
            if (line >> id >> sex >> age >> com ) { //>> empstat) {
                if (com) { // comorbidity default is false, so only need to handle true
                    getPersonByID(id)->setComorbidity(COMORBID);
                }
            }
        }
    }

    if (immunityFilename.length()>0) {
        cerr << "ERROR: Reading in immunity file not currently supported." << endl;
        return false;
/*
        ifstream immiss(immunityFilename);
        if (!immiss) {
            cerr << "ERROR: " << immunityFilename << " not found." << endl;
            return false;
        }
        int part;
        vector<int> parts;
        istringstream line;
        int line_no = 0;
        while ( getline(immiss,buffer) ) {
            line_no++;
            line.clear();
            line.str(buffer);
            while (line >> part) parts.push_back(part);

            // 1+ without age, 2+ with age
            if (parts.size() == 1 + NUM_OF_SEROTYPES or parts.size() == 2 + NUM_OF_SEROTYPES) {
                const int id = parts[0];
                Person* person = getPersonByID(id);
                unsigned int offset = parts.size() - NUM_OF_SEROTYPES;
                vector<pair<int,Serotype> > infection_history;
                for (unsigned int f=offset; f<offset+NUM_OF_SEROTYPES; f++) {
                    Serotype s = (Serotype) (f - offset);
                    const int infection_time = parts[f];
                    if (infection_time == 0) {
                        continue; // no infection for this serotype
                    } else if (infection_time<0) {
                        infection_history.push_back(make_pair(infection_time, s));
                    } else {
                        cerr << "ERROR: Found positive-valued infection time in population immunity file:\n\t";
                        cerr << "person " << person->getID() << ", serotype " << s+1 << ", time " << infection_time << "\n\n";
                        cerr << "Infection time should be provided as a negative integer indicated how many days\n";
                        cerr << "before the start of simulation the infection began.";
                        exit(-359);
                    }
                }
                sort(infection_history.begin(), infection_history.end());
                for (auto p: infection_history) person->infect(p.second, p.first + _day);
            } else if (parts.size() == 0) {
                continue; // skipping blank line, or line that doesn't start with ints
            } else {
                cerr << "ERROR: Unexpected number of values on one line in population immunity file.\n\t";
                cerr << "line num, line: " << line_no << ", " << buffer << "\n\n";
                cerr << "Expected " << 1+NUM_OF_SEROTYPES << " values (person id followed by infection time for each serotype),\n";
                cerr << "found " << parts.size() << endl;
                exit(-361);
            }
            parts.clear();
        }
        immiss.close();
        */
    }

    // keep track of all age cohorts for aging and mortality
    _personAgeCohort.clear();
    _personAgeCohort.resize(NUM_AGE_CLASSES, vector<Person*>(0));

    for (Person* p: _people) {
        int age = p->getAge();
        assert(age<NUM_AGE_CLASSES);
        _personAgeCohort[age].push_back(p);
    }

/*    if (swapFilename == "") {
        _uniformSwap = true;
    } else {
        iss.open(swapFilename);
        if (!iss) {
            cerr << "ERROR: " << swapFilename << " not found." << endl;
            return false;
        }

        int id1, id2;
        double prob;
        istringstream line;

        while ( getline(iss, buffer) ) {
            line.clear();
            line.str(buffer);

            if (line >> id1 >> id2 >> prob) {
                Person* person = getPersonByID(id1);
                if (person) person->appendToSwapProbabilities(make_pair(id2, prob));
            }
        }
        iss.close();
        _uniformSwap = false;
    }*/
    return true;
}


bool Community::loadLocations(string locationFilename,string networkFilename) {
    ifstream iss(locationFilename);
    if (!iss) {
        cerr << "ERROR: locations file " << locationFilename << " not found." << endl;
        return false;
    }
    _location.clear();

    char buffer[500];
    int locID, hfid;
    string locTypeStr;
    string essentialStr;
    double locX, locY;
    istringstream line(buffer);
    map<Location*, int> house_hospitalID_lookup;
    map<int, Location*> hospitalPtr_lookup;

    while (iss) {
        iss.getline(buffer,500);
        line.clear();
        line.str(buffer);
        //if (line >> locID >> locTypeStr >> locX >> locY) {
        //if (line >> locID >> locX >> locY >> locTypeStr >> essentialStr) {
        if (line >> locID >> locX >> locY >> locTypeStr >> essentialStr >> hfid) {
            if (locID != (signed) _location.size()) {
                cerr << "ERROR: Location ID's must be sequential integers starting at 0" << endl;
                cerr << "locID vs _location.size(): " << locID << " " << _location.size() << endl;
                return false;
            }
            const LocationType locType = (locTypeStr == "h") ? HOUSE :
                                         (locTypeStr == "w") ? WORK :
                                         (locTypeStr == "s") ? SCHOOL :
                                         (locTypeStr == "hf") ? HOSPITAL :
                                         (locTypeStr == "n") ? NURSINGHOME :
                                             NUM_OF_LOCATION_TYPES;

            if (locType == NUM_OF_LOCATION_TYPES) {
                cerr << "ERROR: Parsed unknown location type: " << locTypeStr << " from location file: " << locationFilename << endl;
                return false;
            }

            const int essential = (essentialStr == "y" or essentialStr == "NA") ? 1 :
                                  (essentialStr == "n") ? 0 : -1;

            if (essential == -1) {
                cerr << "ERROR: Unknown value for \"essential\" status: " << essentialStr << " from location file: " << locationFilename << endl;
                return false;
            }

            Location* newLoc = new Location();
            newLoc->setX(locX);
            newLoc->setY(locY);
            newLoc->setType(locType); // may be redundant--would save 1 mb per million locations to omit, so probably not worth removing
            if (locType == HOUSE) newLoc->setRiskAversion(gsl_rng_uniform(RNG));
            newLoc->setEssential((bool) essential);
            _location.push_back(newLoc);
            _location_map[locType].insert(newLoc);

            // temp look-up structures so we can quickly map houses to their associated hospitals.
            // those hospital pointers don't necessarily exist until we're done processing locations
            if (locType == HOSPITAL) { hospitalPtr_lookup[newLoc->getID()] = newLoc; }
            if (hfid >= 0) { house_hospitalID_lookup[newLoc] = hfid; } // true for houses, nursinghomes, maybe others in the future
        }
    }
    iss.close();

    for (const auto& kv: house_hospitalID_lookup) { kv.first->setHospital(hospitalPtr_lookup.at(kv.second)); }

    iss.open(networkFilename);
    if (!iss) {
        cerr << "ERROR: network file " << networkFilename << " not found." << endl;
        return false;
    }
    int locID1, locID2;
    while (iss) {
        iss.getline(buffer,500);
        line.clear();
        line.str(buffer);
        if (line >> locID1 >> locID2) { // data (non-header) line
            //      cerr << locID1 << " , " << locID2 << endl;
            assert(locID1 >= 0 and locID2 >= 0);
            assert(locID1 < (signed) _location.size() and locID2 < (signed) _location.size());
            _location[locID1]->addNeighbor(_location[locID2]);            // should check for ID
            _location[locID2]->addNeighbor(_location[locID1]);
        }
    }
    iss.close();

    return true;
}


Person* Community::getPersonByID(int id) {
    // This assumes that IDs start at 1, and tries to guess
    // that person with ID id is in position id-1
    // TODO - make that not true (about starting at 1)
    if(id < 0 or id > (signed) getNumPeople()) {
        cerr << "ERROR: failed to find person with id " << id << " max: " << getNumPeople() << endl;
        assert(id > 0 and id <= (signed) getNumPeople());
    }

    assert (_people[id]->getID() == id);
    return _people[id];
}


// infect - infects person id
Infection* Community::infect(int id) {
    Person* person = getPersonByID(id);
    return person->infect(-1, _date, 0);
}


void Community::vaccinate(CatchupVaccinationEvent cve) {
    // This approach to vaccination is somewhat problematic.  Age classes can be vaccinated multiple times,
    // so the probability of an individual being vaccinated becomes 1 - (1 - ve.coverage)^n, where n is the number
    // of times an age class is specified, either explicitly or implicitly by using a negative value for age

    // Valid coverage and age?
    assert(cve.coverage >= 0.0 and cve.coverage <= 1.0);
    assert(cve.age <= (signed) _personAgeCohort.size());

    for (Person* p: _personAgeCohort[cve.age]) {
        assert(p != NULL);
        if (!p->isVaccinated()
            and cve.coverage > gsl_rng_uniform(RNG)
            and p->isSeroEligible(_par->vaccineSeroConstraint, _par->seroTestFalsePos, _par->seroTestFalseNeg)
           ) {
            p->vaccinate(cve.simDay);
            if (_par->vaccineBoosting or p->getNumVaccinations() < _par->numVaccineDoses) _revaccinate_set.insert(p);
        }
    }
}


void Community::updateVaccination() {
    for (Person* p: _revaccinate_set) {
        if (not p->isVaccinated()) {
            // may be in set unnecessarily because of vaccination before last birthday
            _revaccinate_set.erase(p);
            continue;
        }
        const int timeSinceLastVaccination = p->daysSinceVaccination(_day);
        // TODO: Since updateVaccination only gets called on birthdays, the following only has an effect
        // when the intervals are a multiple of years
        if (p->getNumVaccinations() < _par->numVaccineDoses and timeSinceLastVaccination >= _par->vaccineDoseInterval) {
            // multi-dose vaccination
            p->vaccinate(_day);
            if (p->getNumVaccinations() == _par->numVaccineDoses) _revaccinate_set.erase(p); // we're done
        } else if (_par->vaccineBoosting and timeSinceLastVaccination >= _par->vaccineBoostingInterval) {
            // booster dose
            p->vaccinate(_day);
        }
    }
}


void Community::targetVaccination(Person* p) {
    if (_day < _par->vaccineTargetStartDate) return; // not starting yet
    // expected to be run on p's birthday
    if (p->getAge()==_par->vaccineTargetAge
        and not p->isVaccinated()
        and p->isSeroEligible(_par->vaccineSeroConstraint, _par->seroTestFalsePos, _par->seroTestFalseNeg)
       ) {
        // standard vaccination of target age; vaccinate w/ probability = coverage
        if (gsl_rng_uniform(RNG) < _par->vaccineTargetCoverage) p->vaccinate(_day);
        if (_par->vaccineBoosting or _par->numVaccineDoses > 1) _revaccinate_set.insert(p);
    }
}


vector<pair<size_t,double>> Community::getMeanNumSecondaryInfections() const {
    // this is not how many secondary infections occurred on day=index,
    // but how many secondary infections will ultimately be caused by each person
    // who got infected on day=index
    vector<vector<double>> daily_secondary_infections;
    for (Person* p: _people) {
        for (const Infection* inf: p->getInfectionHistory()) {
            const int infection_onset = inf->getInfectedTime();
            assert(infection_onset >= 0);
            if ((unsigned) infection_onset >= daily_secondary_infections.size()) { daily_secondary_infections.resize(infection_onset+1); }
            // number of secondary infections resulting from the infection that started on this date
            daily_secondary_infections[infection_onset].push_back(inf->secondary_infection_tally());
        }
    }
    vector<pair<size_t, double>> daily_Rt(daily_secondary_infections.size());
    for (size_t day = 0; day < daily_secondary_infections.size(); ++day) {
        daily_Rt[day] = make_pair(daily_secondary_infections[day].size(), mean(daily_secondary_infections[day]));
        cerr << "day, incidence, Rt: " << day << " " << daily_Rt[day].first << " " << daily_Rt[day].second << endl;
    }
    return daily_Rt;
}


void Community::reportCase(int onsetDate, long int reportDate) { // long int b/c reportDate can be a bit greater than max int
    assert(onsetDate >= 0);
    assert(reportDate >= 0);
    // onset == sample collection date; FL doesn't report when symptoms began
    if ((unsigned) onsetDate < _numDetectedCasesOnset.size()) _numDetectedCasesOnset[onsetDate]++;
    if ((unsigned) reportDate < _numDetectedCasesReport.size()) _numDetectedCasesReport[reportDate]++;
}


void Community::reportDeath(int /*eventDate*/, long int reportDate) {
    assert(reportDate >= 0);
    if ((unsigned) reportDate < _numDetectedDeaths.size()) _numDetectedDeaths[reportDate]++;
}


void Community::updatePersonStatus() {
    // TODO - add support for all disease outcomes
    for (Person* p: _people) {
        Location* day_loc = p->getDayLoc();
        if (p->inHospital(_day)) {
            if (p->getHospitalizedTime()==_day) {
                // health care employees may already be at the facility where they would receive treatment
                if (day_loc != p->getHospital()) { p->goToHospital(); }
            }
        } else if (_day > 0 and p->inHospital(_day - 1)) { // they were in hospital yeserday, but no longer
            if (day_loc != p->getHospital()) { p->leaveHospital(); }
        }

        if (p->isSurveilledPerson()) {
            if (p->getNumNaturalInfections() == 0) continue;              // no infection/outcomes to tally
            if (p->getInfectedTime()==_day) _numNewlyInfected[_day]++;
            if (p->getSymptomTime()==_day) {                              // started showing symptoms today
                _numNewlySymptomatic[_day]++;
                if (p->isVaccinated()) {
                    _numVaccinatedCases[_day]++;
                }
            }

            if (p->isSevere(_day)) {
                _numSeverePrev[_day]++;
            }

            if (p->inHospital(_day)) {
                _numHospPrev[_day]++;
                if (p->getHospitalizedTime()==_day) {
                    _numHospInc[_day]++;
                }
                if (p->inIcu(_day)) {
                    _numIcuPrev[_day]++;
                    if (p->getIcuTime()==_day) {
                        _numIcuInc[_day]++;
                    }
                }
            }

            if (p->isNewlyDead(_day)) {
                _numNewlyDead[_day]++;
            }
            /*if (p->getWithdrawnTime()==_day) {                            // started withdrawing
                p->getLocation(HOME_MORNING)->addPerson(p,WORK_DAY);       // stays at home at mid-day
                p->getLocation(WORK_DAY)->removePerson(p,WORK_DAY);        // does not go to work
            } else if (p->isWithdrawn(_day-1) and
            p->getRecoveryTime()==_day) {                                 // just stopped withdrawing
                p->getLocation(WORK_DAY)->addPerson(p,WORK_DAY);           // goes back to work
                p->getLocation(HOME_MORNING)->removePerson(p,WORK_DAY);    // stops staying at home
            }*/
        }
    }
    return;
}


void Community::flagInfectedLocation(LocationType locType, Location* _pLoc, int day) {
    // pass in location type in case we want to model e.g. nursing home workers not going to work, but residents staying put
    assert(day >= 0);
    if ((unsigned) day < _par->runLength) _isHot[day][locType][_pLoc]++;
}


Infection* Community::trace_contact(int &infectee_id, Location* source_loc, int infectious_count) {
    // Identify who was the source of an exposure event (tracing backward)
    // This function currently assumes all co-located infected people are equally likely
    // to be the cause of transmission.  May be something we want to relax in the future.
    const bool loc_is_hospital = source_loc->getType() == HOSPITAL;
    vector<Person*> infected_candidates;
    for (Person* p: source_loc->getPeople()) {
        if(p->isInfectious(_day)
          // possibilities:
          // 1.) We're not looking at a hospital, and this person isn't in a hospital (and thus is here)
          // 2.) We are looking at a hospital, and this person is infected and working at this hospital
          // 3.) We are looking at a hospital, and this person normally works here, but has been admitted here
          // 4.) We are looking at a hospital, and this person normally works here, but has been admitted elsewhere <-- the tricky one
          // 5.) We are looking at a hospital, and this person normally works elsewhere, but has been admitted here
          and (not p->inHospital(_day) or (loc_is_hospital and source_loc == p->getHospital()))
          and not p->isDead(_day)) {
            infected_candidates.push_back(p);
        }
    }

    if ((signed) infected_candidates.size() != infectious_count) {
        cerr << "found vs expected, day " << _day << ": " << infected_candidates.size() << " " << infectious_count << endl;
        cerr << "new lookup value (should == expected): " << _isHot[_day][source_loc->getType()][source_loc] << endl;
        cerr << "Problematic location:\n";
        source_loc->dumper();
        for (Person* p: source_loc->getPeople()) {
            if(p->isInfectious(_day)
              and (not p->inHospital(_day) or (loc_is_hospital and source_loc == p->getHospital()))
              and not p->isDead(_day)) {
                cerr << "\nInfected person " << p->getID() << " normal day loc: " << (p->getDayLoc() ? p->getDayLoc()->getID() : -1) << endl;
                cerr << "Hospitalized?: " << p->inHospital(_day) << endl;
                p->getInfection()->dumper();
            }
        }
    }

    assert((signed) infected_candidates.size() == infectious_count);
    Person* infectee = choice(RNG, infected_candidates);
    infectee_id = infectee->getID();
    return infectee->getInfection();
}


void Community::within_household_transmission() {
    for (const auto hot: _isHot[_day][HOUSE]) {
        Location* loc = hot.first;
        int infectious_count = hot.second;
        //cerr << "\t\t\t\thousehold, count: " << loc->getID() << ", " << infectious_count << endl;
        if (infectious_count > 0) {
            const double T = 1.0 - pow(1.0 - _par->household_transmissibility, infectious_count);
            _transmission(loc, loc->getPeople(), T, infectious_count);
        }
    }
    return;
}


double Community::social_distancing(int _day) {
    return timedInterventions[SOCIAL_DISTANCING][_day];
}


void Community::between_household_transmission() {
    for (auto hot : _isHot[_day][HOUSE]) {
        Location* loc = hot.first;
        const int infectious_count = hot.second;
        // TODO -- it would make more sense to call this risk loving (but that would be weird in this context), instead of risk aversion.
        // this is not intuitive right now
        if (loc->getRiskAversion() > social_distancing(_day)) { // this household is not cautious enough to avoid interactions
            const int hh_size = loc->getNumPeople();
            const float hh_prev = (float) infectious_count / hh_size;
            if (hh_prev > 0.0) {
                for (Location* neighbor: loc->getNeighbors()) {
                    if (neighbor->getRiskAversion() > social_distancing(_day)) {
                        const double T = _par->social_transmissibility * hh_prev;
                        _transmission(loc, neighbor->getPeople(), T, infectious_count);
                    }
                }
            }
        }
    }
    return;
}


void Community::school_transmission() {
    // Transmission for school employees is considered school transmission, not workplace transmission
    for (auto hot : _isHot[_day][SCHOOL]) {
        Location* loc = hot.first;
        const int infectious_count = hot.second;
        const int school_size = loc->getNumPeople();
        if (infectious_count > 0 and school_size > 1) {
            const double T = (1.0 - timedInterventions[SCHOOL_CLOSURE][_day]) * _par->school_transmissibility * infectious_count/(school_size - 1.0);
            _transmission(loc, loc->getPeople(), T, infectious_count);
        }
    }
    return;
}


void Community::workplace_transmission() {
    // Transmission for school employees is considered school transmission, not workplace transmission
    for (auto hot : _isHot[_day][WORK]) {
        Location* loc = hot.first;
        const int infectious_count = hot.second;
        // if non-essential businesses are closed, skip this workplace
        if (loc->isNonEssential() and timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE][_day]) {
            continue;
        }
        const int workplace_size = loc->getNumPeople();
        if (infectious_count > 0 and workplace_size > 1) {
            const double T = (1.0 - social_distancing(_day)) * _par->workplace_transmissibility * infectious_count/(workplace_size - 1.0);
            _transmission(loc, loc->getPeople(), T, infectious_count);
        }
    }
    return;
}


void Community::hospital_transmission() {
    for (auto hot : _isHot[_day][HOSPITAL]) {
        Location* loc = hot.first;
        const int infectious_count = hot.second;
        const int hospital_census = loc->getNumPeople(); // workers + patients
        if (infectious_count > 0 and hospital_census > 1) {
            const double T = _par->hospital_transmissibility * infectious_count/(hospital_census - 1.0);
            _transmission(loc, loc->getPeople(), T, infectious_count);
        }
    }
    return;
}


void Community::nursinghome_transmission() {
    for (auto hot : _isHot[_day][NURSINGHOME]) {
        Location* loc = hot.first;
        const int infectious_count = hot.second;
        const int nursinghome_census = loc->getNumPeople(); // workers + residents
        if (infectious_count > 0 and nursinghome_census > 1) {
            const double T = _par->nursinghome_transmissibility * infectious_count/(nursinghome_census - 1.0);
            _transmission(loc, loc->getPeople(), T, infectious_count);
        }
    }
    return;
}


void Community::_transmission(Location* source_loc, vector<Person*> at_risk_group, const double T, const int infectious_count) {
    for (Person* p: at_risk_group) {
        if (gsl_rng_uniform(RNG) < T) {
            int infectee_id = INT_MIN;
            Infection* infection = nullptr;
            if (_par->traceContacts) {
                infection = trace_contact(infectee_id, source_loc, infectious_count);
            }
            Infection* transmission = p->infect(infectee_id, _date, source_loc->getID()); // infect() tests for whether person is infectable
            if (infection and transmission) { infection->log_transmission(transmission); } // are we contact tracing, and did transmission occur?
        }
    }
}


void Community::updateHotLocations() {
    for (size_t locType = 0; locType < NUM_OF_LOCATION_TYPES; ++locType) {
        _isHot[_day][(LocationType) locType].clear();
    }
}


void Community::tick() {
    _day = _date->day();
    within_household_transmission();
    nursinghome_transmission();
    workplace_transmission();
    if (not timedInterventions[SCHOOL_CLOSURE][_day]) school_transmission();
    between_household_transmission();
// TODO - nursing home interactions

//    if (isWeekday(dow)) {
//        school_transmission();
//    }
//    local_transmission();

    //updateVaccination();

    updatePersonStatus();                                            // make people stay home or return to work
    // Hospital transmission should happen after updating person status,
    // b/c that's when hospitals find out they are receiving a patient.
    // People cannot transmit in a hospital and outside of a hospital on the same day,
    // as transmission for other locations checks whether the person will be admitted
    // on this day.
    hospital_transmission();
    updateHotLocations();

    return;
}


// getNumInfected - counts number of infected residents
size_t Community::getNumInfected(int day) {
    size_t count=0;
    //for (Person* p: _people) { if (p->isInfected(day)) count++; }
    for (Person* p: _people) { count += p->isInfected(day); }
    return count;
}


size_t Community::getNumInfectious(int day) {
    size_t count=0;
    for (Person* p: _people) { if (p->isInfectious(day)) count++; }
    return count;
}


// getNumSymptomatic - counts number of symptomatic residents
size_t Community::getNumSymptomatic(int day) {
    size_t count=0;
    for (Person* p: _people) { if (p->isSymptomatic(day)) count++; }
    return count;
}

// getNumSusceptible - counts number of susceptible residents
size_t Community::getNumNaive() {
    size_t count = 0;
    for (Person* p: _people) {
        if (p->isNaive()) count++;
    }
    return count;
}
