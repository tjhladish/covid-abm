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
vector< map<LocationType, map<Location*, map<double, vector<Person*>>, Location::LocPtrComp>>> Community::_isHot;
set<Person*> Community::_revaccinate_set;
vector<size_t> Community::_numDetectedCasesOnset;
vector<size_t> Community::_numDetectedCasesReport;
vector<size_t> Community::_numDetectedHospitalizations;
vector<size_t> Community::_numDetectedDeaths;
vector<size_t> Community::_cumulIncByOutcome(NUM_OF_OUTCOME_TYPES, 0);

//vector<Person*> Community::_peopleByAge;
//map<int, set<pair<Person*,Person*> > > Community::_delayedBirthdays;

int mod(int k, int n) { return ((k %= n) < 0) ? k+n : k; } // correct for non-negative n


Community::Community(const Parameters* parameters, Date* date) :
//    _exposedQueue(MAX_INCUBATION, vector<Person*>(0)),
    _numNewlyInfected(parameters->runLength), // +1 not needed; runLength is already a valid size
    _numNewVocInfections(parameters->runLength),
    _numNewlySymptomatic(parameters->runLength),
    _numNewlySevere(parameters->runLength),
    _numNewlyCritical(parameters->runLength),
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
    _numDetectedHospitalizations.resize(_par->runLength);
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
    _numNewVocInfections.resize(_par->runLength);
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
    _numNewVocInfections.clear();
    _numNewlySymptomatic.clear();
    _numNewlyDead.clear();
    _numVaccinatedCases.clear();
}


bool Community::loadPopulation(string populationFilename, string comorbidityFilename, string publicActivityFilename, string immunityFilename) {
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
    int pid, hid, age, sex, did;//, empstat;
    while ( getline(iss, buffer) ) {
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
        if (line >> pid >> hid >> sex >> age >> did ) { //>> empstat) {
            if (pid != (signed) _people.size()) { // ensures indexing stays consistent
                cerr << "ERROR: Person ID's must be sequential integers starting at 0" << endl;
                return false;
            }
            assert((signed) _location.size() > hid);
            assert((signed) _location.size() > did);

            Person* p = new Person();
            _people.push_back(p);
            p->setAge(age);
            p->setSex((SexType) sex);
            p->setDaysImmune(_par->sampleDaysImmune(RNG));

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

    if (comorbidityFilename.length() > 0) {
        iss.open(comorbidityFilename);

        if (!iss) {
            cerr << "ERROR: " << comorbidityFilename << " not found." << endl;
            return false;
        }

        bool com;
        while ( getline(iss, buffer) ) {
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
            if (line >> pid >> sex >> age >> com ) { //>> empstat) {
                if (com) { // comorbidity default is false, so only need to handle true
                    getPersonByID(pid)->setComorbidity(COMORBID);
                }
            }
        }
    }

    if (publicActivityFilename.length() > 0) {
        iss.open(publicActivityFilename);

        if (!iss) {
            cerr << "ERROR: " << publicActivityFilename << " not found." << endl;
            return false;
        }

        string buffer;
        istringstream line;

        int pid, locid;
        while ( getline(iss, buffer) ) {
            line.clear();
            line.str(buffer);
            /*
            pid dest_locid_1 dest_locid_2 dest_locid_3 dest_locid_4 dest_locid_5
            0 11321 12073 15279 10654 23786
            1 15042 9939 10607 14346 20246
            2 18748 24527 20246 16895 21665
            */
            if (line >> pid) {
                Person* p = getPersonByID(pid);
                // underage people and LTCF residents do not engage in commercial activities
                if (p->getAge() < 18 or p->getLongTermCare()) { continue; }

                while (line >> locid) {
                    assert((signed) _location.size() > locid);
                    Location* loc = _location[locid];
                    const PublicTransmissionType pub_risk = loc->getPublicTransmissionRisk();
                    assert(pub_risk == LOW_PUBLIC_TRANSMISSION or pub_risk == HIGH_PUBLIC_TRANSMISSION);
                    p->addPatronizedLocation(loc);
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


bool Community::loadLocations(string locationFilename, string networkFilename) {
    ifstream iss(locationFilename);
    if (!iss) {
        cerr << "ERROR: locations file " << locationFilename << " not found." << endl;
        return false;
    }
    _location.clear();

    char buffer[500];
    int locid, hfid;
    string locTypeStr;
    string essentialStr;
    double locX, locY;
    double compliance;
    string publicTransmissionRiskStr;

    istringstream line(buffer);
    map<Location*, int> house_hospitalID_lookup;
    map<int, Location*> hospitalPtr_lookup;

    while (iss) {
        iss.getline(buffer,500);
        line.clear();
        line.str(buffer);
        //if (line >> locid >> locTypeStr >> locX >> locY) {
        //if (line >> locid >> locX >> locY >> locTypeStr >> essentialStr) {
        if (line >> locid >> locX >> locY >> locTypeStr >> essentialStr >> hfid) {
            if (locid != (signed) _location.size()) {
                cerr << "ERROR: Location ID's must be sequential integers starting at 0" << endl;
                cerr << "locid vs _location.size(): " << locid << " " << _location.size() << endl;
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
            newLoc->setEssential((bool) essential);

            if ((line >> compliance) and compliance >= 0) { // predetermined compliance values
                assert(compliance <= 1.0);
                newLoc->setRiskiness(1.0 - compliance);
            } else if (locType == HOUSE) {                  // compliance values not specified, so determined at runtime
                newLoc->setRiskiness(gsl_rng_uniform(RNG));
            }

            PublicTransmissionType public_transmision_risk = NO_PUBLIC_TRANSMISSION;

            if ((line >> publicTransmissionRiskStr) and publicTransmissionRiskStr != "N") {
                public_transmision_risk = (publicTransmissionRiskStr == "H") ? HIGH_PUBLIC_TRANSMISSION :
                                          (publicTransmissionRiskStr == "L") ? LOW_PUBLIC_TRANSMISSION :
                                              NUM_OF_PUBLIC_TRANSMISSION_TYPES;

                if (public_transmision_risk == NUM_OF_PUBLIC_TRANSMISSION_TYPES) {
                    cerr << "ERROR: Unknown value for \"public transmission\" status: " << publicTransmissionRiskStr << " from location file: " << locationFilename << endl;
                    return false;
                } else if (public_transmision_risk != NO_PUBLIC_TRANSMISSION) {
                    _public_locations.push_back(newLoc);
                }
            }

            newLoc->setPublicTransmissionRisk(public_transmision_risk);

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
    int locid1, locid2;
    while (iss) {
        iss.getline(buffer,500);
        line.clear();
        line.str(buffer);
        if (line >> locid1 >> locid2) { // data (non-header) line
            //      cerr << locid1 << " , " << locid2 << endl;
            assert(locid1 >= 0 and locid2 >= 0);
            assert(locid1 < (signed) _location.size() and locid2 < (signed) _location.size());
            _location[locid1]->addNeighbor(_location[locid2]);            // should check for ID
            _location[locid2]->addNeighbor(_location[locid1]);
        }
    }
    iss.close();

    return true;
}


Person* Community::getPersonByID(int pid) {
    // This assumes that IDs start at 1, and tries to guess
    // that person with ID id is in position id-1
    // TODO - make that not true (about starting at 1)
    if(pid < 0 or pid > (signed) getNumPeople()) {
        cerr << "ERROR: failed to find person with id " << pid << " max: " << getNumPeople() << endl;
        assert(pid > 0 and pid <= (signed) getNumPeople());
    }

    assert (_people[pid]->getID() == pid);
    return _people[pid];
}


// infect - infects person id
Infection* Community::infect(int id, StrainType strain) {
    Person* person = getPersonByID(id);
    return person->infect(_date, strain);
}


void Community::vaccinate(CatchupVaccinationEvent cve) {
    // This approach to vaccination is somewhat problematic.  Age classes can be vaccinated multiple times,
    // so the probability of an individual being vaccinated becomes 1 - (1 - ve.coverage)^n, where n is the number
    // of times an age class is specified, either explicitly or implicitly by using a negative value for age

    // Valid coverage and age?
    assert(cve.coverage >= 0.0 and cve.coverage <= 1.0);
    assert(cve.age <= _personAgeCohort.size());

    const size_t duration = cve.campaignDuration;
    for (Person* p: _personAgeCohort[cve.age]) {
        assert(p != NULL);
        // is today this person's vaccination day (during a multi-day campaign)?
        if (p->getID() % duration == _day % duration) {
            if (!p->isVaccinated()
                and cve.coverage > gsl_rng_uniform(RNG)
                and p->isSeroEligible(_par->vaccineSeroConstraint, _par->seroTestFalsePos, _par->seroTestFalseNeg)
               ) {
                const bool success = p->vaccinate(_day);
                if (success and (_par->vaccineBoosting or (p->getNumVaccinations() < _par->numVaccineDoses))) {
                //cerr << _par->vaccineBoosting << " | " << p->getNumVaccinations() << " < " << _par->numVaccineDoses << endl;
                _revaccinate_set.insert(p); }
            }
        }
    }
}


void Community::updateVaccination() {
    vector<Person*> to_erase;
    for (Person* p: _revaccinate_set) {
        if (not p->isVaccinated()) {
            // may be in set unnecessarily because of vaccination before last birthday
            to_erase.push_back(p);
            continue;
        }
        const int timeSinceLastVaccination = p->daysSinceVaccination(_day);
        if (p->getNumVaccinations() < _par->numVaccineDoses and timeSinceLastVaccination >= _par->vaccineDoseInterval) {
            // multi-dose vaccination
            p->vaccinate(_day);
            if (p->getNumVaccinations() == _par->numVaccineDoses) { to_erase.push_back(p); } // we're done
        } else if (_par->vaccineBoosting and timeSinceLastVaccination >= _par->vaccineBoostingInterval) {
            // booster dose
            p->vaccinate(_day);
        }
    }

    for (Person* p: to_erase) { _revaccinate_set.erase(p); }
}


void Community::targetVaccination(Person* p) {
    if (_day < _par->vaccineTargetStartDate) return; // not starting yet
    // expected to be run on p's birthday
    if (p->getAge()==_par->vaccineTargetAge
        and not p->isVaccinated()
        and p->isSeroEligible(_par->vaccineSeroConstraint, _par->seroTestFalsePos, _par->seroTestFalseNeg)
       ) {
        // standard vaccination of target age; vaccinate w/ probability = coverage
        if (gsl_rng_uniform(RNG) < _par->vaccineTargetCoverage) { p->vaccinate(_day); }
        if (_par->vaccineBoosting or _par->numVaccineDoses > 1) { _revaccinate_set.insert(p); }
    }
}


vector<pair<size_t,double>> Community::getMeanNumSecondaryInfections() const {
    // this is not how many secondary infections occurred on day=index,
    // but how many secondary infections will ultimately be caused by each person
    // who got infected on day=index
    vector<vector<double>> daily_secondary_infections(_par->runLength);
    for (Person* p: _people) {
        for (const Infection* inf: p->getInfectionHistory()) {
            const int infection_onset = inf->getInfectedTime();
            if (infection_onset < 0) { continue; } // historical infection
            // number of secondary infections resulting from the infection that started on this date
            daily_secondary_infections[infection_onset].push_back(inf->secondary_infection_tally());
            // Uncomment this to do offspring distribution/dispersion analyses
            //cerr << "secondary: " << infection_onset << " " << inf->secondary_infection_tally() << endl;
        }
    }
    vector<pair<size_t, double>> daily_Rt(daily_secondary_infections.size(), {0, 0.0});
    for (size_t day = 0; day < daily_secondary_infections.size(); ++day) {
        if (daily_secondary_infections[day].size()) {
            daily_Rt[day] = make_pair(daily_secondary_infections[day].size(), mean(daily_secondary_infections[day]));
        }
//        cerr << "day, incidence, Rt: " << day << " " << daily_Rt[day].first << " " << daily_Rt[day].second << endl;
    }
    return daily_Rt;
}


void Community::reportCase(int onsetDate, long int reportDate, bool hospitalized) { // long int b/c reportDate can be a bit greater than max int
    assert(onsetDate >= 0);
    assert(reportDate >= 0);
    // onset == sample collection date; FL doesn't report when symptoms began
    if ((unsigned) onsetDate < _numDetectedCasesOnset.size()) { _numDetectedCasesOnset[onsetDate]++; }
    if ((unsigned) reportDate < _numDetectedCasesReport.size()) {
        _numDetectedCasesReport[reportDate]++;
        // it's not clear exactly how to interpret the date on which the state reports a hospitalization
        if (hospitalized) { _numDetectedHospitalizations[reportDate]++; }
    }
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
            if (p->getNumNaturalInfections() == 0) {
                continue;              // no infection/outcomes to tally
            } else {
                // the methods used with Infection below generally are available for Person, but this should be faster
                const Infection* inf = p->getInfection();
                if (inf->getInfectedTime()==_day) {
                    _numNewlyInfected[_day]++;
                    if (inf->getStrain() != WILDTYPE) { _numNewVocInfections[_day]++; }
                }

                if (inf->getSymptomTime()==_day) {                              // started showing symptoms today
                    _numNewlySymptomatic[_day]++;
                    if (p->isVaccinated()) { _numVaccinatedCases[_day]++; }
                }

                if (inf->isSevere(_day)) { _numSeverePrev[_day]++; }

                if (inf->inHospital(_day)) {
                    _numHospPrev[_day]++;
                    if (inf->getHospitalizedTime()==_day) {
                        _numHospInc[_day]++;
                    }
                    if (inf->inIcu(_day)) {
                        _numIcuPrev[_day]++;
                        if (inf->getIcuTime()==_day) {
                            _numIcuInc[_day]++;
                        }
                    }
                }

                if (inf->getSevereTime()==_day)   { _numNewlySevere[_day]++; }
                if (inf->getCriticalTime()==_day) { _numNewlyCritical[_day]++; }
                if (p->isNewlyDead(_day))       { _numNewlyDead[_day]++; }
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
    }
    return;
}


void Community::flagInfectedLocation(Person* person, double relInfectiousness, LocationType locType, Location* _pLoc, int day) {
    assert(day >= 0);
    if ((unsigned) day < _par->runLength) _isHot[day][locType][_pLoc][relInfectiousness].push_back(person);
}


Infection* Community::trace_contact(Person* &infecter, Location* source_loc, const map<double, vector<Person*>> &infectious_groups) {
    // Identify who was the source of an exposure event (tracing backward)
    // First we determine which group did the infecting (grouped by infectiousness),
    // then we chose the person within the group who is the infecter
    vector<double> individual_weights;
    vector<double> group_weights;
    double total = 0.0;
    for (const auto& [relInfectiousness, people]: infectious_groups) {
        individual_weights.push_back(relInfectiousness);
        const double group_weight = relInfectiousness * people.size();
        total += group_weight;
        group_weights.push_back(group_weight);
    }

    double r = total * gsl_rng_uniform(RNG);
    size_t idx;
    for (idx = 0; idx<group_weights.size(); ++idx) {
        if (r < group_weights[idx]) {
            break;
        } else {
            r -= group_weights[idx];
        }
    }
    infecter = choice(RNG, infectious_groups.at(individual_weights[idx]));

    // sanity check to make sure we've found a legit candidate
    const vector<Person*> people = source_loc->getPeople();
    assert(infecter->isInfectious(_day)
            and (not infecter->inHospital(_day) or (source_loc->getType() == HOSPITAL and source_loc == infecter->getHospital()))
            and not infecter->isDead(_day)
            and find(people.begin(), people.end(), infecter) != people.end());

    return infecter->getInfection();
}


double Community::social_distancing(int _day) {
    return timedInterventions[SOCIAL_DISTANCING][_day];
}


double _tally_infectiousness (const map<double, vector<Person*>> infectious_groups) {
    double infectious_weight = 0.0;
    for (const auto& [relInfectiousness, people]: infectious_groups) {
        infectious_weight += relInfectiousness * people.size();
    }
    return infectious_weight;
}


void Community::within_household_transmission() {
    for (const auto& [loc, infectious_groups]: _isHot[_day][HOUSE]) {
        const double infectious_weight = _tally_infectiousness(infectious_groups);
        const double hazard =  _par->household_transmissibility * _par->seasonality_on(_date) * infectious_weight;
        const double T = 1.0 - exp(-hazard);
        _transmission(loc, loc->getPeople(), infectious_groups, T);
    }
    return;
}


void Community::between_household_transmission() {
    for (const auto& [loc, infectious_groups]: _isHot[_day][HOUSE]) {
        const double infectious_weight = _tally_infectiousness(infectious_groups);
        // ↓↓↓ this model made it almost impossible to stop transmission using SD
        //if (social_distancing(_day) - loc->getRiskiness() < gsl_rng_uniform(RNG)) { // this household is not cautious enough to avoid interactions

        // if people are riskier than the current SD level, they interact with friends
        // if they are less risky than current SD, they may do so, depending on how much more cautious they are
        // ↓↓↓ this household is not cautious enough to avoid interactions
        if (loc->getRiskiness() > social_distancing(_day)) {
            const int hh_size = loc->getNumPeople();
            for (Location* neighbor: loc->getNeighbors()) {
                // ↓↓↓ this line needs to match the model above, with neighbor in for loc
                if (neighbor->getRiskiness() > social_distancing(_day)) {
                    const double hazard = _par->social_transmissibility * _par->seasonality_on(_date) * infectious_weight / hh_size;
                    const double T = 1.0 - exp(-hazard);
                    _transmission(loc, neighbor->getPeople(), infectious_groups, T);
                }
            }
        }
    }
    return;
}


void Community::workplace_transmission() {
    // Transmission for school employees is considered school transmission, not workplace transmission
    // This includes all other employees, as well as consumer visits to restaurants, bars, retail locations, and religious facilities
    for (const auto& [loc, infectious_groups]: _isHot[_day][WORK]) {
        // if non-essential businesses are closed, skip this workplace
        const int workplace_size = loc->getNumPeople() + loc->getNumVisitors();
        if (workplace_size < 2 or (loc->isNonEssential() and timedInterventions[NONESSENTIAL_BUSINESS_CLOSURE][_day])) {
            continue;
        }

        const double infectious_weight = _tally_infectiousness(infectious_groups);
        const PublicTransmissionType pt_risk = loc->getPublicTransmissionRisk();

        if (infectious_weight > 0) {
            const double hazard = _par->workplace_transmissibility
                                  // TODO -- see if we can find a way to motivate how much more risky high risk places are
                                  // TODO -- check to see if high risk places actually are causing 4x as much transmission as other workplaces
                                  //* (pt_risk == HIGH_PUBLIC_TRANSMISSION ? 4.0 : (1.0 - social_distancing(_day))*0.25) // 4.0 and 0.25 b/c/ of 80/20 rule
                                  * (pt_risk == HIGH_PUBLIC_TRANSMISSION ? 4.0 : 0.25) // 4.0 and 0.25 b/c/ of 80/20 rule
                                  * _par->seasonality_on(_date)
                                  * infectious_weight/(workplace_size - 1.0);

            const double T = 1.0 - exp(-hazard);
            vector<Person*> all_people = loc->getVisitors();
            const vector<Person*> workers = loc->getPeople();
            all_people.insert( all_people.end(), workers.begin(), workers.end() );
            _transmission(loc, all_people, infectious_groups, T);
        }
    }
    return;
}


void Community::school_transmission() {
    // Transmission for school employees is considered school transmission, not workplace transmission
    const double hazard_coef = (1.0 - timedInterventions[SCHOOL_CLOSURE][_day]) * _par->school_transmissibility * _par->seasonality_on(_date);
    if (hazard_coef != 0.0) {
        for (const auto& [loc, infectious_groups]: _isHot[_day][SCHOOL]) {
            const int school_size = loc->getNumPeople();
            if (school_size < 2) { continue; }
            const double infectious_weight = _tally_infectiousness(infectious_groups);
            const double hazard = hazard_coef * infectious_weight/(school_size - 1.0);
            const double T = 1.0 - exp(-hazard);
            _transmission(loc, loc->getPeople(), infectious_groups, T);
        }
    }
    return;
}


void Community::hospital_transmission() {
    for (const auto& [loc, infectious_groups]: _isHot[_day][HOSPITAL]) {
        const int hospital_census = loc->getNumPeople(); // workers + patients
        if (hospital_census < 2) { continue; }
        const double infectious_weight = _tally_infectiousness(infectious_groups);
        const double hazard = _par->hospital_transmissibility * _par->seasonality_on(_date) * infectious_weight/(hospital_census - 1.0);
        const double T = 1.0 - exp(-hazard);
        _transmission(loc, loc->getPeople(), infectious_groups, T);
    }
    return;
}

/*
generic_location_transmission(_isHot[_day][HOSPITAL], _par->hospital_transmissibility * _par->seasonality_on(_date));
TODO - switch to this?
void Community::generic_location_transmission(const auto& hot_location_type_data, const double base_T) {
    for (const auto& [loc, infectious_groups]: hot_location_type_data) {
        const int census = loc->getNumPeople(); // workers + patients
        if (census < 2) { continue; }
        const double infectious_weight = _tally_infectiousness(infectious_groups);
        const double hazard = base_T * infectious_weight/(census - 1.0);
        const double T = 1.0 - exp(-hazard);
        _transmission(loc, loc->getPeople(), infectious_groups, T);
    }
    return;
}*/


void Community::nursinghome_transmission() {
    for (const auto& [loc, infectious_groups]: _isHot[_day][NURSINGHOME]) {
        const int nursinghome_census = loc->getNumPeople(); // workers + residents
        if (nursinghome_census < 2) { continue; }
        const double infectious_weight = _tally_infectiousness(infectious_groups);
        const double hazard = _par->nursinghome_transmissibility * _par->seasonality_on(_date) * infectious_weight/(nursinghome_census - 1.0);
        const double T = 1.0 - exp(-hazard);
        _transmission(loc, loc->getPeople(), infectious_groups, T);
    }
    return;
}


void Community::_transmission(Location* source_loc, vector<Person*> at_risk_group, const map<double, vector<Person*>> &infectious_groups, const double T) {
    const bool check_susceptibility = false;
    for (Person* p: at_risk_group) {
        if (gsl_rng_uniform(RNG) < T and p->isInfectable(_date->day())) { // isInfectable() seems likely to be more expensive
            Person* infecter     = nullptr;
            Infection* source_infection = nullptr;
            // because we now support multiple strains, we always have to trace, in order to determine what the infecting strain would be
            source_infection = trace_contact(infecter, source_loc, infectious_groups);
            Infection* transmission = p->infect(infecter, _date, source_loc, source_infection->getStrain(), check_susceptibility); // infect() tests for whether person is infectable
            if (source_infection and transmission) { source_infection->log_transmission(transmission); } // did we contact trace, and did transmission occur?
        }
    }
}


void Community::updateHotLocations() {
    for (size_t locType = 0; locType < NUM_OF_LOCATION_TYPES; ++locType) {
        _isHot[_day][(LocationType) locType].clear();
    }
}


void Community::public_activity() {
    for (Person* p: _people) {
        const vector<Location*> locs = p->getPatronizedLocations();
        const int num_locs = locs.size();
        const bool avoiding_high_risk_places = p->getRiskiness() < social_distancing(_day);
        if (num_locs and not p->inHospital(_day)) { // does this person patronize businesses?
            const float lambda_visit_hours_per_day = 1.0;
            size_t visit_hours = gsl_ran_poisson(RNG, lambda_visit_hours_per_day);
            map<Location*, int> visit_map;

            for (size_t i = 0; i < visit_hours; ++i) {
                size_t idx = gsl_rng_uniform_int(RNG, num_locs); // sampled on [0, num_locs-1]
                const Location* loc = locs[idx];
                const PublicTransmissionType risk_lvl = loc->getPublicTransmissionRisk();
                // if this person is being cautious, they avoid high risk environments
                if (avoiding_high_risk_places and (risk_lvl == HIGH_PUBLIC_TRANSMISSION)) {
                    continue;
                } else {
                    visit_map[locs[idx]]++; // and if they aren't, they don't
                }
            }

            for (const auto& [loc, hours]: visit_map) {
                loc->addVisitor(p, hours);

                if (p->isInfectious(_day)) {
                    double relInfectiousness = p->getRelInfectiousness();
                    relInfectiousness *= (double) hours/8; // infectiousness compared to an employee with an 8h shift
                    flagInfectedLocation(p, relInfectiousness, loc->getType(), loc, _day);
                }
            }
        }
    }
}


void Community::clear_public_activity() {
    for (Location* loc: _public_locations) {
        loc->clearVisitors();
    }
}


void Community::tick() {
    _day = _date->day();
    within_household_transmission();
    between_household_transmission();

    public_activity();
    workplace_transmission();
    clear_public_activity();

    if (not timedInterventions[SCHOOL_CLOSURE][_day]) school_transmission();
    nursinghome_transmission();

//    if (isWeekday(dow)) {
//        school_transmission();
//    }
//    local_transmission();

    updateVaccination();

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
