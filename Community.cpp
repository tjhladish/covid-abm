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

const Parameters* Community::_par;
vector< map<LocationType, map<Location*, int, Location::LocPtrComp>>> Community::_isHot;
set<Person*> Community::_revaccinate_set;
//vector<Person*> Community::_peopleByAge;
//map<int, set<pair<Person*,Person*> > > Community::_delayedBirthdays;

int mod(int k, int n) { return ((k %= n) < 0) ? k+n : k; } // correct for non-negative n

Community::Community(const Parameters* parameters) :
//    _exposedQueue(MAX_INCUBATION, vector<Person*>(0)),
    _numNewlyInfected(parameters->runLength), // +1 not needed; runLength is already a valid size
    _numNewlySymptomatic(parameters->runLength),
    _numVaccinatedCases(parameters->runLength),
    _numSevereCases(parameters->runLength)
    {
    _par = parameters;
    _day = 0;
    //_mortality = NULL;
    //_bNoSecondaryTransmission = false;
    //_uniformSwap = true;
//    for (int a = 0; a<NUM_AGE_CLASSES; a++) _personAgeCohortSizes[a] = 0;
    _isHot.resize(_par->runLength);
    for (auto &e: _isHot) {
        for (size_t locType = 0; locType < NUM_OF_LOCATION_TYPES; ++locType) {
            e[(LocationType) locType] = {};
        }
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
    _numVaccinatedCases.clear();

//    _exposedQueue.resize(MAX_INCUBATION, vector<Person*>(0));
    _numNewlyInfected.resize(_par->runLength);
    _numNewlySymptomatic.resize(_par->runLength);
    _numVaccinatedCases.resize(_par->runLength);
}


Community::~Community() {
    if (_people.size() > 0) { for (Person* p: _people) delete p; }

    Person::reset_ID_counter();

    for (auto &e: _isHot) e.clear();

    for (unsigned int i = 0; i < _location.size(); i++ ) delete _location[i];
    _location.clear();

    for (auto& kv: _location_map) {
        kv.second.clear();
    }
    _location_map.clear();

//    _exposedQueue.clear();
    _personAgeCohort.clear();
    _numNewlyInfected.clear();
    _numNewlySymptomatic.clear();
    _numVaccinatedCases.clear();
}


bool Community::loadPopulation(string populationFilename, string immunityFilename) {
    ifstream iss(populationFilename.c_str());

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
//            p->setLocation(_location[hid], HOME);
//            p->setLocation(_location[did], DAY);
            _location[hid]->addPerson(p); // TODO -- currently locations that can be a home *AND* a workplace cannot exist
            if (did >= 0) {
                p->setDayLoc(_location[did]); // currently just any non-home daytime location--may be a school; -1 == NA
                _location[did]->addPerson(p); // distinguish between the reasons why someone is there without querying the person
            }
            if (_location[hid]->getType() == NURSINGHOME) p->setLongTermCare(true);
            //assert(age<NUM_AGE_CLASSES);
 //           agecounts[age]++;
        }
    }
    iss.close();

    //_peopleByAge = _people;
    //sort(_peopleByAge.begin(), _peopleByAge.end(), PerPtrComp());

    if (immunityFilename.length()>0) {
        cerr << "ERROR: Reading in immunity file not currently supported." << endl;
        return false;
/*
        ifstream immiss(immunityFilename.c_str());
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
//    _personAgeCohort.clear();
//    _personAgeCohort.resize(NUM_AGE_CLASSES, vector<Person*>(0));

/*    for (Person* p: _people) {
        int age = p->getAge();
        assert(age<NUM_AGE_CLASSES);
        _personAgeCohort[age].push_back(p);
        _personAgeCohortSizes[age]++;
    }*/

/*    if (swapFilename == "") {
        _uniformSwap = true;
    } else {
        iss.open(swapFilename.c_str());
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
    ifstream iss(locationFilename.c_str());
    if (!iss) {
        cerr << "ERROR: locations file " << locationFilename << " not found." << endl;
        return false;
    }
    _location.clear();

    char buffer[500];
    int locID;
    string locTypeStr;
    string essentialStr;
    double locX, locY;
    istringstream line(buffer);

    while (iss) {
        iss.getline(buffer,500);
        line.clear();
        line.str(buffer);
        //if (line >> locID >> locTypeStr >> locX >> locY) {
        if (line >> locID >> locX >> locY >> locTypeStr >> essentialStr) {
            if (locID != (signed) _location.size()) {
                cerr << "ERROR: Location ID's must be sequential integers starting at 0" << endl;
                return false;
            }
            const LocationType locType = (locTypeStr == "h") ? HOUSE :
                                         (locTypeStr == "w") ? WORK :
                                         (locTypeStr == "s") ? SCHOOL :
                                         (locTypeStr == "t") ? HOSPITAL :
                                         (locTypeStr == "n") ? NURSINGHOME :
                                             NUM_OF_LOCATION_TYPES;

            if (locType == NUM_OF_LOCATION_TYPES) {
                cerr << "ERROR: Parsed unknown location type: " << locTypeStr << " from location file: " << locationFilename << endl;
                return false;
            }

            const int essential = (essentialStr == "y" or essentialStr == "NA") ? 1 :
                                  (essentialStr == "n") ? 0 :
                                      -1;
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
        }
    }
    iss.close();

    iss.open(networkFilename.c_str());
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
bool Community::infect(int id, int day) {
    Person* person = getPersonByID(id);

    bool result =  person->infect(-1, day, 0);
    if (result) _numNewlyInfected[_day]++;
    return result;
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


void Community::updateDiseaseStatus() {
    // TODO - add support for all disease outcomes
    for (Person* p: _people) {
        if (p->getNumNaturalInfections() == 0) continue;
        if (p->getSymptomTime()==_day) {                              // started showing symptoms today
            _numNewlySymptomatic[_day]++;
            if (p->isVaccinated()) {
                _numVaccinatedCases[_day]++;
            }
            if (p->isSevere(_day)) {                          // symptoms will be severe at onset
                _numSevereCases[_day]++;     // if they're going to be severe
            }
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
    return;
}


void Community::flagInfectedLocation(LocationType locType, Location* _pLoc, int day) {
    // pass in location type in case we want to model e.g. nursing home workers not going to work, but residents staying put
    //if (day < _par->runLength) _isHot[day].insert(_pLoc);
    if (day < _par->runLength) _isHot[day][locType][_pLoc]++;
}


void Community::within_household_transmission() {
    //for (Location* loc: _location_map[HOUSE]) { // TODO -- tracking 'hot' households will avoid looping through the vast majority
    //for (Location* loc: _isHot[_day][HOUSE]) {
    for (const auto hot: _isHot[_day][HOUSE]) {
        Location* loc = hot.first;
        int infectious_count = hot.second;
//        for (Person* p: loc->getPeople()) { // if isHot is a map with a count, we wouldn't need this loop at all
//            infectious_count += p->isInfectious(_day) and not p->inHospital(_day);
//        }
        if (infectious_count > 0) {
            const double T = 1.0 - pow(1.0 - _par->household_transmissibility, infectious_count);
            for (Person* p: loc->getPeople()) {
                if (gsl_rng_uniform(RNG) < T) {
                    p->infect((int) HOUSE, _day, loc->getID()); // infect() tests for whether person is infectable
                }
            }
        }
    }
    return;
}


double Community::social_distancing(int _day) const {
    return _par->timedInterventionEffect(SOCIAL_DISTANCING, _day);
}


void Community::between_household_transmission() {
    for (Location* loc: _location_map[HOUSE]) { // TODO -- tracking 'hot' households will avoid looping through the vast majority
        if (_isHot[_day][HOUSE].count(loc) == 0) continue;
        if (loc->getRiskAversion() > social_distancing(_day)) { // this household is not cautious enough to avoid interactions
            int infectious_count = 0;
            const int hh_size = loc->getNumPeople();

            for (Person* p: loc->getPeople()) { // if isHot is a map with a count, we wouldn't need this loop at all
                infectious_count += p->isInfectious(_day) and not p->isSevere(_day);
            }

            const float hh_prev = (float) infectious_count / hh_size;
            if (hh_prev > 0.0) {
                for (Location* neighbor: loc->getNeighbors()) {
                    if (neighbor->getRiskAversion() > social_distancing(_day)) {
                        const double T = _par->social_transmissibility * hh_prev;
                        for (Person* p: loc->getPeople()) {
                            if (gsl_rng_uniform(RNG) < T) {
                                p->infect((int) HOUSE, _day, loc->getID()); // infect() tests for whether person is infectable
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}


void Community::workplace_and_school_transmission() {
    // some (== non-essential) workplaces may be closed, but schools are either all open or all closed
    //location_transmission(_location_map[WORK]);
    //if (not _par->timedInterventionEffect(SCHOOL_CLOSURE, _day)) location_transmission(_location_map[SCHOOL]);
    location_transmission(_isHot[_day][WORK]);
    if (not _par->timedInterventionEffect(SCHOOL_CLOSURE, _day)) location_transmission(_isHot[_day][SCHOOL]);
}


void Community::location_transmission(map<Location*, int, Location::LocPtrComp> &locations) {
    // TODO -- right now, school and workplace transmission are handled separately
    // Transmission for school employees is considered school transmission, not workplace transmission
    for (auto hot : locations) {
        Location* loc = hot.first;
        const int infectious_count = hot.second;
        // if non-essential businesses are closed, skip this workplace
        if (loc->isNonEssential() and _par->timedInterventionEffect(NONESSENTIAL_BUSINESS_CLOSURE, _day)) {
            continue;
        }
        const int workplace_size = loc->getNumPeople();
        if (infectious_count > 0 and workplace_size > 1) {
            const double T = _par->workplace_transmissibility * infectious_count/(workplace_size - 1.0);
            for (Person* p: loc->getPeople()) {
                if (gsl_rng_uniform(RNG) < T) {
                    p->infect((int) HOME, _day, loc->getID()); // infect() tests for whether person is infectable
                }
            }
        }
    }
    return;
}


void Community::location_transmission(set<Location*, Location::LocPtrComp> &locations) {
    // TODO -- right now, school and workplace transmission are handled separately
    // Transmission for school employees is considered school transmission, not workplace transmission
    for (Location* loc: locations) { // TODO -- track 'hot' workplaces/schools
        // if non-essential businesses are closed, skip this workplace
        if (loc->isNonEssential() and _par->timedInterventionEffect(NONESSENTIAL_BUSINESS_CLOSURE, _day)) {
            continue;
        }

        int infectious_count = 0;
        for (Person* p: loc->getPeople()) { // if isHot is a map with a count, we wouldn't need this loop at all
            infectious_count += p->isInfectious(_day) and not p->isSevere(_day);
        }
        const int workplace_size = loc->getNumPeople();
        if (infectious_count > 0 and workplace_size > 1) {
            const double T = _par->workplace_transmissibility * infectious_count/(workplace_size - 1.0);
            for (Person* p: loc->getPeople()) {
                if (gsl_rng_uniform(RNG) < T) {
                    p->infect((int) HOME, _day, loc->getID()); // infect() tests for whether person is infectable
                }
            }
        }
    }
    return;
}


void Community::updateHotLocations() {
    for (size_t locType = 0; locType < NUM_OF_LOCATION_TYPES; ++locType) {
        _isHot[_day][(LocationType) locType] = {};
    }
}


void Community::tick(int day) {
    _day = day;
    within_household_transmission();
    workplace_and_school_transmission(); // does not currently include schools or churches
    between_household_transmission();
// TODO - nursing home interactions

//    if (isWeekday(dow)) {
//        school_transmission();
//    }
//    local_transmission();

    //updateVaccination();

    updateDiseaseStatus();                                            // make people stay home or return to work
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
