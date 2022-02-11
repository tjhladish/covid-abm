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
#include "Vac_Campaign.h"

using namespace covid::standard;
using covid::util::mean;
using covid::util::uniform_choice;
using covid::util::weighted_choice;
using covid::util::choose_k;
using covid::util::merge_vectors;

const Parameters* Community::_par;

int mod(int k, int n) { return ((k %= n) < 0) ? k+n : k; } // correct for non-negative n

Community::Community(const Parameters* parameters, Date* date) :
    _numNewlyInfected(parameters->runLength, 0), // +1 not needed; runLength is already a valid size
    _numNewlySymptomatic(parameters->runLength, 0),
    _numNewlySevere(parameters->runLength, 0),
    _numNewlyCritical(parameters->runLength, 0),
    _numNewlyDead(parameters->runLength, 0),
    _numVaccinatedCases(parameters->runLength, 0),
    _numSeverePrev(parameters->runLength, 0),
    _numHospInc(parameters->runLength, 0),
    _numHospPrev(parameters->runLength, 0),
    _numIcuInc(parameters->runLength, 0),
    _numIcuPrev(parameters->runLength, 0),
    _numDetectedCasesOnset(parameters->runLength, 0),
    _numDetectedCasesReport(parameters->runLength, 0),
    _numDetectedHospitalizations(parameters->runLength, 0),
    //_numDetectedDeaths(parameters->runLength, 0),
    _numDetectedDeathsOnset(parameters->runLength, 0),
    _numDetectedDeathsReport(parameters->runLength, 0),
    _cumulIncByOutcome(NUM_OF_OUTCOME_TYPES, 0),
    _isHot(parameters->runLength)
    {
    _par = parameters;
    _date = date;
    _day = 0;
    for (int strain = 0; strain < (int) NUM_OF_STRAIN_TYPES; ++strain) {
        _numNewInfectionsByStrain[(StrainType) strain] = vector<size_t>(_par->runLength);
    }
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
    for (int strain = 0; strain < (int) NUM_OF_STRAIN_TYPES; ++strain) {
        _numNewInfectionsByStrain[(StrainType) strain] = vector<size_t>(_par->runLength);
    }
    _numNewlySymptomatic.resize(_par->runLength);
    _numNewlyDead.resize(_par->runLength);
    _numVaccinatedCases.resize(_par->runLength);
}


Community::~Community() {
    if (_date) delete _date;
    if (_people.size() > 0) { for (Person* p: _people) delete p; }
    if (vac_campaign) { delete vac_campaign; }

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
    _numNewInfectionsByStrain.clear();
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
            p->setStartingNaturalEfficacy(_par->sampleStartingNaturalEfficacy(RNG));
            p->setImmunityQuantile(gsl_rng_uniform(RNG));
            p->setNaturalImmunityDuration(_par->immunityDuration(p->getImmunityQuantile(), p->getStartingNaturalEfficacy()));
            //p->setDaysImmune(_par->sampleDaysImmune(RNG));

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

double _calculatePixel(double coord) {
    return (floor(coord / 0.01) * 0.01) + 0.005;
}


bool Community::loadLocations(string locationFilename, string networkFilename) {
    ifstream iss(locationFilename);
    if (!iss) {
        cerr << "ERROR: locations file " << locationFilename << " not found." << endl;
        return false;
    }
    _location.clear();

    string buffer;

    int locid, hfid;
    string locTypeStr;
    string essentialStr;
    double locX, locY;
    double compliance;
    string publicTransmissionRiskStr;

    istringstream line;

    map<Location*, int> house_hospitalID_lookup;
    map<int, Location*> hospitalPtr_lookup;

    while ( getline(iss, buffer) ) {
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

            double xPixel = _calculatePixel(locX);
            double yPixel = _calculatePixel(locY);
            newLoc->setPixel(xPixel, yPixel);
            _pixelMap[{xPixel, yPixel}].push_back(newLoc);

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

    while ( getline(iss, buffer) ) {
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
    return person->infect(this, _date, strain);
}


void Community::vaccinate() {
    Vaccinee* v = vac_campaign->next_vaccinee(_day);
    while (v) { // v is not nullptr, e.g. there is actually a dose and a person who might be vaccinatable
        if (((!v->get_person()->isVaccinated() and v->get_person()->isSeroEligible()) // unvaccinated & eligible
          or (v->get_status() == REVACCINATE_QUEUE)) // or scheduled for a subsequent dose
          and vac_campaign->vaccinate(v, _day)) { // and person isn't dead, so got vaccinated
            vac_campaign->tally_dose(_day, v); // tally dose used, and person vaccinated

            // multi-dose vaccines
            if (v->getNumVaccinations() < _par->numVaccineDoses) {
                vac_campaign->schedule_revaccination(_day + _par->vaccineDoseInterval, v);
            }

            // vaccines that require regular boosting
            if (_par->vaccineBoosting) {
                vac_campaign->schedule_revaccination(_day + _par->vaccineBoostingInterval, v);
            }
        }
        delete v;
        v = vac_campaign->next_vaccinee(_day);
    }

    // if we ran out of doses for today but still have people that need revaccination, reschedule them for tomorrow
    vac_campaign->reschedule_remaining_revaccinations(_day);
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


//void Community::reportDeath(int /*eventDate*/, long int reportDate) {
//    assert(reportDate >= 0);
//    if ((unsigned) reportDate < _numDetectedDeaths.size()) _numDetectedDeaths[reportDate]++;
//}
void Community::reportDeath(int onsetDate, long int reportDate) {
    assert(onsetDate >= 0);
    assert(reportDate >= 0);
    if ((unsigned) onsetDate < _numDetectedDeathsOnset.size()) { _numDetectedDeathsOnset[onsetDate]++; }
    if ((unsigned) reportDate < _numDetectedDeathsReport.size()) { _numDetectedDeathsReport[reportDate]++; }
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
                    _numNewInfectionsByStrain.at(inf->getStrain())[_day]++;
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
        if (p->isQuarantining(_day) and not p->isQuarantining(_day+1)) {
            p->endQuarantine();
        }
    }
    return;
}


void Community::flagInfectedLocation(Person* person, double relInfectiousness, LocationType locType, Location* _pLoc, int day) {
    assert(day >= 0);
    if ((unsigned) day < _par->runLength) _isHot[day][locType][_pLoc][relInfectiousness].push_back(person);
}


//vector<double> trans_type(3, 0.0); // asymptomatic, presymptomatic, symptomatic, for logging transmission type
Infection* Community::trace_contact(Person* &infecter, Location* source_loc, const map<double, vector<Person*>> &infectious_groups) {
    // Identify who was the source of an exposure event (tracing backward)
    // First we determine which group did the infecting (grouped by infectiousness),
    // then we choose the person within the group who is the infecter
    vector<double> relInfectiousnessValues;
    vector<double> group_weights;
    double total = 0.0;
    for (const auto& [relInfectiousness, people]: infectious_groups) {
        relInfectiousnessValues.push_back(relInfectiousness);
        const double group_weight = relInfectiousness * people.size();
        total += group_weight;
        group_weights.push_back(group_weight);
    }

    size_t group_idx = weighted_choice(RNG, group_weights);
    infecter = uniform_choice(RNG, infectious_groups.at(relInfectiousnessValues[group_idx]));

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
        const double hazard =  _par->household_transmission_haz_mult * _par->seasonality_on(_date) * infectious_weight;
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
                    const double hazard = _par->social_transmission_haz_mult * _par->seasonality_on(_date) * infectious_weight / hh_size;
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

        const double high_pt_risk_haz_mult = 4.0;
        const double norm_pt_risk_haz_mult = 0.25;

        if (infectious_weight > 0) {
            const double hazard = _par->workplace_transmission_haz_mult
                                  // TODO -- see if we can find a way to motivate how much more risky high risk places are
                                  // TODO -- check to see if high risk places actually are causing 4x as much transmission as other workplaces
                                  //* (pt_risk == HIGH_PUBLIC_TRANSMISSION ? 4.0 : (1.0 - social_distancing(_day))*0.25) // 4.0 and 0.25 b/c/ of 80/20 rule
                                  * (pt_risk == HIGH_PUBLIC_TRANSMISSION ? high_pt_risk_haz_mult : norm_pt_risk_haz_mult) // 4.0 and 0.25 b/c/ of 80/20 rule
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
    const double hazard_coef = (1.0 - timedInterventions[SCHOOL_CLOSURE][_day]) * _par->school_transmission_haz_mult * _par->seasonality_on(_date);
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
        const double hazard = _par->hospital_transmission_haz_mult * _par->seasonality_on(_date) * infectious_weight/(hospital_census - 1.0);
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
        const double hazard = _par->nursinghome_transmission_haz_mult * _par->seasonality_on(_date) * infectious_weight/(nursinghome_census - 1.0);
        const double T = 1.0 - exp(-hazard);
        _transmission(loc, loc->getPeople(), infectious_groups, T);
    }
    return;
}

void Community::_transmission(Location* source_loc, vector<Person*> at_risk_group, const map<double, vector<Person*>> &infectious_groups, const double T) {
    const bool check_susceptibility = true;
    for (Person* p: at_risk_group) {
        if (p->isQuarantining(_date->day())) { continue; }
        if (gsl_rng_uniform(RNG) < T) {
            Person* infecter     = nullptr;
            Infection* source_infection = nullptr;
            // because we now support multiple strains, we always have to trace, in order to determine what the infecting strain would be
            source_infection = trace_contact(infecter, source_loc, infectious_groups);
            // infect() tests for whether person is infectable
            Infection* transmission = p->infect(this, infecter, _date, source_loc, source_infection->getStrain(), check_susceptibility);
            if (source_infection and transmission) {
                source_infection->log_transmission(transmission);
                // for logging transmission type
                // if (infecter->isSymptomatic(_day)) {
                //     trans_type[2]++;
                // } else if (infecter->getInfection()->symptomatic()) {
                //     trans_type[1]++;
                // } else {
                //     trans_type[0]++;
                // }
           } // did we contact trace, and did transmission occur?
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

void _conditional_insert(set<Person*> &into, const vector<Person*> &from, set<Person*> &ref, const size_t day) {
    for(Person* p : from) {
        if(p->isAlive(day) and ref.insert(p).second) { into.insert(p); }
    }
}

vector< set<Person*> > Community::traceForwardContacts() {
    //vector< set<Person*> > tracedContacts(_par->contactTracingDepth);      // returned data structure of traced contacts of the priumary cases categorized by depth
    vector< set<Person*> > tracedContacts(_par->contactTracingDepth + 1);      // returned data structure of traced contacts of the priumary cases categorized by depth
    // only do contact tracing after the start date set in main
    if(_day < _par->beginContactTracing) { return tracedContacts; }

    set<Person*> tracedCases;
    for(Person* p : _people) {
        if(p->isAlive(_day) and p->hasBeenInfected()) {
            Infection* mostRecentInfection = p->getInfection();
            if(mostRecentInfection->isDetectedOn(_day) and gsl_rng_uniform(RNG) < _par->contactTracingCoverage){
                tracedCases.insert(p);
            }
        }
    }

    tracedContacts[0] = tracedCases;

    set<Person*> allTracedPeople;               // used to ensure all people are traced only once
    set<Person*> allInterviewedPeople;          // used to ensure all people are interviewed only once

    //for(size_t depth = 0; depth < _par->contactTracingDepth; ++depth) {
    for(size_t depth = 0; depth < _par->contactTracingDepth; ++depth) {
        //set<Person*> peopleToInterview = depth == 0 ? tracedCases : tracedContacts[depth-1];
        set<Person*> peopleToInterview = depth == 0 ? tracedCases : tracedContacts[depth];

        for(Person* p : peopleToInterview) {
            // set.insert().second returns bool: T if inserted, F if not inserted
            // if TRUE --> person has not been traced before
            if(p->isAlive(_day) and allInterviewedPeople.insert(p).second) {
                // if home is a nursing home, select contacts from other residents + employees
                // number of contacts to find based on poisson distribution with defined expected value
                // all residents + employees of this nursing home
                // find nursing home contacts for this interviewed person
                if(p->getHomeLoc()->getType() == NURSINGHOME) {
                    const size_t numNursingHomeContacts  = gsl_ran_poisson(RNG, _par->contactTracingEV[NURSINGHOME]);
                    vector<Person*> nursingHomeResidents = p->getHomeLoc()->getPeople();
                    vector<Person*> nursingHomeContacts  = (numNursingHomeContacts <= nursingHomeResidents.size())
                                                               ? choose_k(RNG, nursingHomeResidents, numNursingHomeContacts)
                                                               : nursingHomeResidents;

                    // insert contacts into final data structure at this depth
                    //_conditional_insert(tracedContacts[depth], nursingHomeContacts, allTracedPeople, _day);
                    _conditional_insert(tracedContacts[depth + 1], nursingHomeContacts, allTracedPeople, _day);
                } else {
                    // all home residents are known contacts of this interviewed person
                    vector<Person*> allHomeContacts = p->getHomeLoc()->getPeople();
                    //_conditional_insert(tracedContacts[depth], allHomeContacts, allTracedPeople, _day);
                    _conditional_insert(tracedContacts[depth + 1], allHomeContacts, allTracedPeople, _day);

                    // contact trace neighbors
                    const size_t numNeighborContacts = gsl_ran_poisson(RNG, _par->contactTracingEV[HOME]);
                    vector<Person*> allNeighbors;
                    for(Location* loc : p->getHomeLoc()->getNeighbors()) {
                        for(Person* n : loc->getPeople()) { allNeighbors.push_back(n); }
                    }
                    vector<Person*> neighborContacts = (numNeighborContacts <= allNeighbors.size())
                                                           ? choose_k(RNG, allNeighbors, numNeighborContacts)
                                                           : allNeighbors;
                    //_conditional_insert(tracedContacts[depth], neighborContacts, allTracedPeople, _day);
                    _conditional_insert(tracedContacts[depth + 1], neighborContacts, allTracedPeople, _day);

                    // contact trace day location people (work or school)
                    vector<Person*> dayLocContacts;
                    if(p->getDayLoc()) {
                        LocationType lt = p->getDayLoc()->getType();
                        const size_t numDayLocContacts = gsl_ran_poisson(RNG, _par->contactTracingEV[lt]);
                        vector<Person*> dayLocPeople = p->getDayLoc()->getPeople();
                        dayLocContacts = (numDayLocContacts <= dayLocPeople.size()) ? choose_k(RNG, dayLocPeople, numDayLocContacts)
                                                                                    : dayLocPeople;
                    }
                    //_conditional_insert(tracedContacts[depth], dayLocContacts, allTracedPeople, _day);
                    _conditional_insert(tracedContacts[depth + 1], dayLocContacts, allTracedPeople, _day);
                }
            }
        }
// size_t num_quarantined = 0;
//        for (Person* p : tracedContacts[depth]) {
//            if ((not p->isQuarantining(_day)) and (gsl_rng_uniform(RNG) < _par->quarantineProbability[depth])) {
//                p->selfQuarantine(_day, _par->selfQuarantineDuration);
//                num_quarantined++;
//            }
//        }
// cerr << "DEBUG NUM CONTACTS AT DEPTH         " << depth << " IS " << tracedContacts[depth].size() << endl;
// cerr << "DEBUG PROB OF QUARANTINING AT DEPTH " << depth << " IS " << _par->quarantineProbability[depth] << endl;
// cerr << "DEBUG NUM QUARANTINED FROM DEPTH    " << depth << " IS " << num_quarantined << endl;
    }
 size_t num_quarantined = 0;
 for (size_t depth = 0; depth < tracedContacts.size(); ++depth) {
        for (Person* p : tracedContacts[depth]) {
            if ((not p->isQuarantining(_day)) and (gsl_rng_uniform(RNG) < _par->quarantineProbability[depth])) {
                p->selfQuarantine(_day, _par->selfQuarantineDuration);
                num_quarantined++;
            }
        }
 //cerr << "DEBUG NUM CONTACTS AT DEPTH         " << depth << " IS " << tracedContacts[depth].size() << endl;
 //cerr << "DEBUG PROB OF QUARANTINING AT DEPTH " << depth << " IS " << _par->quarantineProbability[depth] << endl;
 //cerr << "DEBUG NUM QUARANTINED FROM DEPTH    " << depth << " IS " << num_quarantined << endl;
 }
    return tracedContacts;
}

void Community::tick() {
    _day = _date->day();

//cerr << endl << "D (u s)    " << vac_campaign->get_doses_available(_day, URGENT_ALLOCATION) << ' ' << vac_campaign->get_doses_available(_day, STANDARD_ALLOCATION) << endl;
//cerr <<         "Q (u s re) " << vac_campaign->get_urgent_queue_size() << ' ' << vac_campaign->get_standard_queue_size() << ' ' << vac_campaign->get_revaccinate_queue_size(_day) << endl;
    if (vac_campaign) { vaccinate(); }
//cerr <<         "V (u s re) " << vac_campaign->get_dose_tally(_day, URGENT_QUEUE) << ' ' << vac_campaign->get_dose_tally(_day, STANDARD_QUEUE) << ' ' <<  vac_campaign->get_dose_tally(_day, REVACCINATE_QUEUE) << endl;

    within_household_transmission();
    between_household_transmission();

    public_activity(); // must come before workplace transmission
    workplace_transmission();
    clear_public_activity();

    if (_date->isWeekday() and not timedInterventions[SCHOOL_CLOSURE][_day]) {
        school_transmission();
    }
    nursinghome_transmission();

    //updateVaccination();

    updatePersonStatus();                                            // make people stay home or return to work
    // Hospital transmission should happen after updating person status,
    // b/c that's when hospitals find out they are receiving a patient.
    // People cannot transmit in a hospital and outside of a hospital on the same day,
    // as transmission for other locations checks whether the person will be admitted
    // on this day.
    hospital_transmission();
    updateHotLocations();

    // do contact tracing to a given depth using reported cases from today if _day is at or after the start of contact tracing
    vector< set<Person*> > tracedContactsByDepth = traceForwardContacts();

    if(vac_campaign) { vac_campaign->reactive_strategy(_day, tracedContactsByDepth, this); } // if there is no reactive strategy, nothing happens

    // output transmission type data
    //if (_day == (int) _par->runLength -1) {
    //    double all_trans = accumulate(trans_type.begin(), trans_type.end(), 0.0);
    //    cerr << "never symptomatic, pre-symptomatic, symptomatic transmission: " << trans_type[0]/all_trans << ", " << trans_type[1]/all_trans << ", " << trans_type[2]/all_trans << endl;
    //}
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

// getNumNaive - counts number of residents with no infection/vaccination history whatsoever
size_t Community::getNumNaive() {
    size_t count = 0;
    for (Person* p: _people) { if (p->isImmuneState(NAIVE)) count++; }
    return count;
}

double Community::doSerosurvey(const ImmuneStateType ist, vector<Person*> &pop, int time) {
    if (pop.size() == 0) { pop = _people; }
    double seropos = 0;
    double seroneg = 0;

    for (Person* p : pop) {
        switch (ist) {
            case NATURAL: // similar(ish) to an N IgG assay
                if (p->hasBeenInfected()) {
                    const Infection* last_inf = p->getInfectionHistory().back();
                    vector<int> possible_last_infection_end_dates = {last_inf->getInfectiousEndTime(), last_inf->getSymptomEndTime()};
                    const int last_infection_end_date = covid::util::max_element(possible_last_infection_end_dates);
                    const int time_since_last_infection = time - last_infection_end_date;
                    double remaining_natural_efficacy = _par->remainingEfficacy(p->getStartingNaturalEfficacy(), time_since_last_infection);

                    if (remaining_natural_efficacy > _par->seroPositivityThreshold) { ++seropos; }
                    else { ++seroneg; }
                } else {
                    ++seroneg;
                }
                break;
            case VACCINATED: // similar(ish) to an S IgG & N IgG assays, (positive and negative, respectively)
                if (p->isVaccinated()) {

                } else {

                }
                break;
            case NATURAL_AND_VACCINATED: // similar(ish) to an S IgG assay
                if (p->hasBeenInfected() or p->isVaccinated()) {

                } else {

                }
                break;
            default:
                break;
        }
    }

    return seropos / (seropos + seroneg);
}

/*
Each month, this function will be called to conduct a hSAR survey for 2 months prior (to increase the chance that we capture all infections an index causes in their household at a point in time)
    Each person infected during the survey month + if no one else in the house was infected within 10 days prior, will be included as an index
    Each index will be queried to see how many secondary infections they caused over the next 2 months within their household
*/
double Community::getHouseholdSecondaryAttackRate(std::vector<Person*> &pop) {
    if (pop.size() == 0) { pop = _people; }     // if no specific population is provided, by default, use the entire population
    if (not (_date->month() >= 3)) { return 0.0; }     // can only begin to survey starting on the 3rd month of the simulation

    int current_sim_day = _day;                         // will be the first day of the month (ensured by calling scope)

    size_t current_julian_month = _date->julianMonth();
    size_t survey_julain_month    = (((int) current_julian_month - 2) >= 1) ? current_julian_month - 2 : ((int) current_julian_month - 2) + 12;
    size_t inbetween_julian_month = (((int) current_julian_month - 1) >= 1) ? current_julian_month - 1 : ((int) current_julian_month - 1) + 12;

    size_t days_in_survey_month    = (_date->isLeap()) ? LEAP_DAYS_IN_MONTH.at(survey_julain_month - 1) : COMMON_DAYS_IN_MONTH.at(survey_julain_month - 1);
    size_t days_in_inbetween_month = (_date->isLeap()) ? LEAP_DAYS_IN_MONTH.at(inbetween_julian_month - 1) : COMMON_DAYS_IN_MONTH.at(inbetween_julian_month - 1);

    assert(current_sim_day - (int) days_in_inbetween_month - (int) days_in_survey_month >= 0);
    int start_simday_survey_month = current_sim_day - (int) days_in_inbetween_month - (int) days_in_survey_month;
    int end_simday_survey_month   = current_sim_day + days_in_survey_month;

    vector<double> individual_SAR_measurements;

    for (Person* p : pop) {
        int secondary_household_infections = 0;
        bool include_in_survey = false;

        if (p->hasBeenInfected() and (p->getInfectedTime() >= start_simday_survey_month and p->getInfectedTime() <= end_simday_survey_month)) {
            for (Person* fam : p->getHomeLoc()->getPeople()) {
                if (fam == p) { continue; }
                if (fam->hasBeenInfected() and (fam->getInfectedTime() >= p->getInfectedTime() - 10 and fam->getInfectedTime() <= p->getInfectedTime())) {
                    include_in_survey = false;
                } else {
                    include_in_survey = true;
                }
            }

            if (not include_in_survey) { continue; }

            Infection* inf = p->getInfection();     // most recent infection
            for (Infection* secondary_inf : inf->get_infections_caused()) {
                if (secondary_inf->getInfectionOwner()->getHomeLoc() == p->getHomeLoc()) { ++secondary_household_infections; }
            }
        }
        if (secondary_household_infections > 0) {
            double SAR = (double) secondary_household_infections / p->getHomeLoc()->getNumPeople();
            individual_SAR_measurements.push_back(SAR);
        }
    }
    return covid::util::mean(individual_SAR_measurements);
}

map<string, double> Community::calculate_daily_direct_VE() {
    map<string, double> VE_map;
    if (vac_campaign) {
        size_t num_ppl_fully_vaxd = 0;
        size_t num_breakthru_infs = 0, num_breakthru_dis = 0, num_breakthru_hosp = 0, num_breakthru_dths = 0, num_breakthru_reported = 0;
        size_t num_unvaxd_infs = 0, num_unvaxd_dis  = 0, num_unvaxd_hosp = 0, num_unvaxd_dths = 0, num_unvaxd_reported = 0;

        for (Person* p : _people) {
            bool fully_vaxd_w_protection = p->isVaccinated() and p->getNumVaccinations() > 1 and p->daysSinceVaccination(_day) >= _par->vaccine_dose_to_protection_lag;
            if (fully_vaxd_w_protection) { ++num_ppl_fully_vaxd; }
            if (p->hasBeenInfected()) {
                if (p->getInfectedTime() == _day) {      // for VEs
                    if (fully_vaxd_w_protection) { ++num_breakthru_infs; }
                    else { ++num_unvaxd_infs; }
                }

                if (p->getSymptomTime() == _day) {       // for VEp
                    if (fully_vaxd_w_protection) { ++num_breakthru_dis; }
                    else { ++num_unvaxd_dis; }
                }

                if (p->getHospitalizedTime() == _day) {  // for VEh
                    if (fully_vaxd_w_protection) { ++num_breakthru_hosp; }
                    else { ++num_unvaxd_hosp; }
                }

                if (p->getDeathTime() == _day) {         // for VEd
                    if (fully_vaxd_w_protection) { ++num_breakthru_dths; }
                    else { ++num_unvaxd_dths; }
                }

                if (p->getInfection()->isDetectedOn(_day)) {
                    if (fully_vaxd_w_protection) { ++num_breakthru_reported; }
                    else { ++num_unvaxd_reported; }
                }
            }
        }

        size_t num_ppl_unvaxd = getNumPeople() - num_ppl_fully_vaxd;
        VE_map["coverage"]    = (double) num_ppl_fully_vaxd / getNumPeople();

        double unvax_inf_risk = (double) num_unvaxd_infs/num_ppl_unvaxd;
        double vax_inf_risk   = (double) num_breakthru_infs/num_ppl_fully_vaxd;
        double inf_risk_ratio = (double) vax_inf_risk/unvax_inf_risk;
        VE_map["VES"]         = 1.0 - inf_risk_ratio;

        double unvax_dis_risk = (double) num_unvaxd_dis/num_ppl_unvaxd;
        double vax_dis_risk   = (double) num_breakthru_dis/num_ppl_fully_vaxd;
        double dis_risk_ratio = (double) vax_dis_risk/unvax_dis_risk;
        VE_map["VEP"]         = 1.0 - dis_risk_ratio;

        double unvax_hosp_risk = (double) num_unvaxd_hosp/num_ppl_unvaxd;
        double vax_hosp_risk   = (double) num_breakthru_hosp/num_ppl_fully_vaxd;
        double hosp_risk_ratio = (double) vax_hosp_risk/unvax_hosp_risk;
        VE_map["VEH"]          = 1.0 - hosp_risk_ratio;

        double unvax_dth_risk = (double) num_unvaxd_dths/num_ppl_unvaxd;
        double vax_dth_risk   = (double) num_breakthru_dths/num_ppl_fully_vaxd;
        double dth_risk_ratio = (double) vax_dth_risk/unvax_dth_risk;
        VE_map["VED"]         = 1.0 - dth_risk_ratio;

        //VE_map["breakthruRatio"] = (double) num_breakthru_infs/(num_breakthru_infs + num_unvaxd_infs);
        VE_map["breakthruRatio"] = (double) num_breakthru_reported/(num_breakthru_reported + num_unvaxd_reported);
        VE_map["vaxInfs"]        = (double) num_breakthru_infs;
        VE_map["unvaxInfs"]      = (double) num_unvaxd_infs;
        VE_map["vaxHosp"]        = (double) num_breakthru_hosp;
        VE_map["unvaxHosp"]      = (double) num_unvaxd_hosp;
    }
    return VE_map;
}
