// Community.h
// Manages relationship between people and locations
#ifndef __COMMUNITY_H
#define __COMMUNITY_H
#include <string>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <numeric>
#include <cmath>
#include <algorithm>
#include "Date.h"
#include "Location.h"
#include "Vac_Campaign.h"

class Person;
// class Vac_Campaign;

// We use this to created a vector of people, sorted by decreasing age.  Used for aging/immunity swapping.
//struct PerPtrComp { bool operator()(const Person* A, const Person* B) const { return A->getAge() > B->getAge(); } };

class CommunityLedger {
    friend class Community;
public:
    CommunityLedger(const Parameters* par) :
        _numNewlyInfected(par->runLength, 0), // +1 not needed; runLength is already a valid size
        _numNewlySymptomatic(par->runLength, 0),
        _numNewlySevere(par->runLength, 0),
        _numNewlyCritical(par->runLength, 0),
        _numNewlyDead(par->runLength, 0),
        _numVaccinatedCases(par->runLength, 0),
        _numSeverePrev(par->runLength, 0),
        _numHospInc(par->runLength, 0),
        _numHospPrev(par->runLength, 0),
        _numIcuInc(par->runLength, 0),
        _numIcuPrev(par->runLength, 0),
        _numDetectedCasesOnset(par->runLength, 0),
        _numDetectedCasesReport(par->runLength, 0),
        _numDetectedHospitalizations(par->runLength, 0),
        //_numDetectedDeaths(parameters->runLength, 0),
        _numDetectedDeathsOnset(par->runLength, 0),
        _numDetectedDeathsReport(par->runLength, 0),
        _cumulIncByOutcome(NUM_OF_OUTCOME_TYPES, 0),
        _isHot(par->runLength),
        _timedInterventions(par->timedInterventions)
    {
        for (int strain = 0; strain < (int) NUM_OF_STRAIN_TYPES; ++strain) {
            _numNewInfectionsByStrain[(StrainType) strain] = vector<size_t>(par->runLength);
        }

        vector<string> inf_by_loc_keys = {"home", "social", "work_staff", "patron", "school_staff", "student", "hcw", "patient", "ltcf_staff", "ltcf_resident"};
        for (string key : inf_by_loc_keys) {
            _numNewlyInfectedByLoc[key] = vector<size_t>(par->runLength, 0);
        }

        for (auto &e: _isHot) {
            for (size_t locType = 0; locType < NUM_OF_LOCATION_TYPES; ++locType) {
                e[(LocationType) locType] = {};
            }
        }
    }

    CommunityLedger(const CommunityLedger& o) {
        _numNewlyInfected = o._numNewlyInfected;
        _numNewInfectionsByStrain = o._numNewInfectionsByStrain;
        _numNewlyInfectedByLoc = o._numNewlyInfectedByLoc;
        _numNewlySymptomatic = o._numNewlySymptomatic;
        _numNewlySevere = o._numNewlySevere;
        _numNewlyCritical = o._numNewlyCritical;
        _numNewlyDead = o._numNewlyDead;
        _numVaccinatedCases = o._numVaccinatedCases;
        _numSeverePrev = o._numSeverePrev;
        _numHospInc = o._numHospInc;
        _numHospPrev = o._numHospPrev;
        _numIcuInc = o._numIcuInc;
        _numIcuPrev = o._numIcuPrev;
        _numDetectedCasesOnset = o._numDetectedCasesOnset;
        _numDetectedCasesReport = o._numDetectedCasesReport;
        _numDetectedHospitalizations = o._numDetectedHospitalizations;
        _numDetectedDeathsOnset = o._numDetectedDeathsOnset;
        _numDetectedDeathsReport = o._numDetectedDeathsReport;
        _cumulIncByOutcome = o._cumulIncByOutcome;
        _isHot = o._isHot;
        _timedInterventions = o._timedInterventions;
    }

    ~CommunityLedger() {   // is any of this necessay?
        _numNewlyInfected.clear();
        _numNewInfectionsByStrain.clear();
        _numNewlyInfectedByLoc.clear();
        _numNewlySymptomatic.clear();
        _numNewlySevere.clear();
        _numNewlyCritical.clear();
        _numNewlyDead.clear();
        _numVaccinatedCases.clear();
        _numSeverePrev.clear();
        _numHospInc.clear();
        _numHospPrev.clear();
        _numIcuInc.clear();
        _numIcuPrev.clear();
        _numDetectedCasesOnset.clear();
        _numDetectedCasesReport.clear();
        _numDetectedHospitalizations.clear();
        _numDetectedDeathsOnset.clear();
        _numDetectedDeathsReport.clear();
        _cumulIncByOutcome.clear();
        _isHot.clear();
        _timedInterventions.clear();
    }
private:
    std::vector<size_t> _numNewlyInfected;
    std::map<StrainType, std::vector<size_t>> _numNewInfectionsByStrain;
    std::map<std::string, std::vector<size_t>> _numNewlyInfectedByLoc;
    std::vector<size_t> _numNewlySymptomatic;                              // true cases, no lag due to detection
    std::vector<size_t> _numNewlySevere;                                   // true cases, no lag due to detection
    std::vector<size_t> _numNewlyCritical;                                 // true cases, no lag due to detection
    std::vector<size_t> _numNewlyDead;                                     // true cases, no lag due to detection
    std::vector<size_t> _numVaccinatedCases;
    std::vector<size_t> _numSeverePrev;
    std::vector<size_t> _numHospInc;
    std::vector<size_t> _numHospPrev;
    std::vector<size_t> _numIcuInc;
    std::vector<size_t> _numIcuPrev;
    std::vector<size_t> _numDetectedCasesOnset;
    std::vector<size_t> _numDetectedCasesReport;
    std::vector<size_t> _numDetectedHospitalizations;
    std::vector<size_t> _numDetectedDeathsOnset;
    std::vector<size_t> _numDetectedDeathsReport;
    std::vector<size_t> _cumulIncByOutcome;

    // groups of infectious people: indexed by day, location type, location ptr, and relative infectiousness
    std::vector< std::map<LocationType, std::map<Location*, std::map<double, std::vector<Person*>>, Location::LocPtrComp>>> _isHot;

    std::map<TimedIntervention, std::vector<double>> _timedInterventions;
};

class Community {
    public:
        Community(const Parameters* parameters, Date* date);
        Community(const Community& o) {
            _par                         = o._par;// is static const
            _date                        = o._date ? new Date(*(o._date)) : nullptr;
            _people                      = std::vector<Person*>(o._people.size());
            _personAgeCohort             = o._personAgeCohort;
            _location                    = std::vector<Location*>(o._location.size());
            _public_locations            = o._public_locations;
            _location_map                = o._location_map;
            // _exposedQueue                = o._exposedQueue;
            _day                         = o._day;
            // _numNewlyInfected            = o._numNewlyInfected;
            // _numNewInfectionsByStrain    = o._numNewInfectionsByStrain;
            // _numNewlyInfectedByLoc       = o._numNewlyInfectedByLoc;
            // _numNewlySymptomatic         = o._numNewlySymptomatic;
            // _numNewlySevere              = o._numNewlySevere;
            // _numNewlyCritical            = o._numNewlyCritical;
            // _numNewlyDead                = o._numNewlyDead;
            // _numVaccinatedCases          = o._numVaccinatedCases;
            // _numSeverePrev               = o._numSeverePrev;
            // _numHospInc                  = o._numHospInc;
            // _numHospPrev                 = o._numHospPrev;
            // _numIcuInc                   = o._numIcuInc;
            // _numIcuPrev                  = o._numIcuPrev;
            // _numDetectedCasesOnset       = o._numDetectedCasesOnset;
            // _numDetectedCasesReport      = o._numDetectedCasesReport;
            // _numDetectedHospitalizations = o._numDetectedHospitalizations;
            //_numDetectedDeaths           = o._numDetectedDeaths;
            // _numDetectedDeathsOnset      = o._numDetectedDeathsOnset;
            // _numDetectedDeathsReport     = o._numDetectedDeathsReport;
            // _cumulIncByOutcome           = o._cumulIncByOutcome;
            cmty_ledger                  = o.cmty_ledger ? new CommunityLedger(*(o.cmty_ledger)) : nullptr;
            // _isHot                       = std::vector< std::map<LocationType, std::map<Location*, std::map<double, std::vector<Person*>>, Location::LocPtrComp>>>(o._isHot.size());
            for (auto &e: cmty_ledger->_isHot) {
                for (size_t locType = 0; locType < NUM_OF_LOCATION_TYPES; ++locType) {
                    e[(LocationType) locType] = {};
                }
            }

            // _peopleByAge                 = o._peopleByAge;
            vac_campaign                 = o.vac_campaign ? new Vac_Campaign(*(o.vac_campaign)) : nullptr;
            // _revaccinate_set             = o._revaccinate_set;
            // timedInterventions           = o.timedInterventions;

            // allocate new location and person objects
            map<Location*, Location*> location_ptr_map;
            for (size_t i = 0; i < o._location.size(); ++i) {
                _location[i] = new Location(*(o._location[i]));
                location_ptr_map[o._location[i]] = _location[i];
            }

            map<Person*, Person*> person_ptr_map;
            for (size_t i = 0; i < o._people.size(); ++i) {
                _people[i] = new Person(*(o._people[i]));
                person_ptr_map[o._people[i]] = _people[i];
            }

            // update pointers to new Location and Person objects
            for(Location* &loc : _location) {
                loc->setHospital(location_ptr_map[loc->getHospital()]);

                vector<Person*> people = loc->getPeople();
                for(Person* &p : people) { p = person_ptr_map[p]; }
                loc->setPeople(people);

                vector<Person*> visitors = loc->getVisitors();
                for(Person* &v : visitors) { v = person_ptr_map[v]; }
                loc->setVisitors(visitors);

                set<Location*, Location::LocPtrComp> neighbors;
                for(Location* n : loc->getNeighbors()) { neighbors.insert(location_ptr_map[n]); }
                loc->setNeighbors(neighbors);
            }

            for(Person* &p : _people) {
                p->setHomeLoc(location_ptr_map[p->getHomeLoc()]);
                p->setDayLoc(location_ptr_map[p->getDayLoc()]);

                vector<Location*> patronizedLocations = p->getPatronizedLocations();
                for(Location* &loc : patronizedLocations) { loc = location_ptr_map[loc]; }
                p->setPatronizedLocations(patronizedLocations);

                for (Infection* inf : p->getInfectionHistory()) {
                    inf->setInfectedPlace(location_ptr_map[inf->getInfectedPlace()]);
                    inf->setInfectionOwner(person_ptr_map[inf->getInfectionOwner()]);
                    inf->setInfectedBy(person_ptr_map[inf->getInfectedBy()]);
                }
            }

            for(auto& v : _personAgeCohort) {
                for(Person* &p : v) { p = person_ptr_map[p]; }
            }

            for(Location* &loc : _public_locations) { loc = location_ptr_map[loc]; }

            // for(auto& kv : _location_map) {
            //     for(Location* &loc : kv.second) { loc = location_ptr_map[loc]; }
            // }

            // groups of infectious people: indexed by day, location type, location ptr, and relative infectiousness
            // std::vector< std::map<LocationType, std::map<Location*, std::map<double, std::vector<Person*>>, Location::LocPtrComp>>> _isHot;
            for (size_t day = 0; day < o.cmty_ledger->_isHot.size(); ++day) {
                // std::map<LocationType, std::map<Location*, std::map<double, std::vector<Person*>>, Location::LocPtrComp>>
                for (const auto& [lt, locMap] : o.cmty_ledger->_isHot[day]) {
                    // std::map<Location*, std::map<double, std::vector<Person*>>, Location::LocPtrComp>
                    for (const auto& [loc, relinfMap] : locMap) {
                        // std::map<double, std::vector<Person*>>
                        for (const auto& [relinfness, v] : relinfMap) {
                            // std::vector<Person*>
                            for (Person* p : v) {
                                Location* new_loc = location_ptr_map[loc];
                                Person* new_per   = person_ptr_map[p];
                                cmty_ledger->_isHot[day][lt][new_loc][relinfness].push_back(new_per);
                            }
                        }
                    }
                }
            }

            if (vac_campaign) {
                Vaccinee_Pool psv = vac_campaign->get_potential_vaccinees(STANDARD_QUEUE);
                Vaccinee_Pool puv = vac_campaign->get_potential_vaccinees(URGENT_QUEUE);

                Eligibility_Q other_sq = o.vac_campaign->get_eligibility_queue(STANDARD_QUEUE);
                Eligibility_Q other_uq = o.vac_campaign->get_eligibility_queue(URGENT_QUEUE);

                Eligibility_Q sq = vac_campaign->get_eligibility_queue(STANDARD_QUEUE);
                Eligibility_Q uq = vac_campaign->get_eligibility_queue(URGENT_QUEUE);
                sq.clear(); sq.resize(_par->numVaccineDoses);
                uq.clear(); uq.resize(_par->numVaccineDoses);

                for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                    for (int bin : vac_campaign->get_unique_age_bins()) {
                        for (Person* &p : psv[dose][bin]) { p = person_ptr_map[p]; }
                        for (Person* &p : puv[dose][bin]) { p = person_ptr_map[p]; }
                    }

                    while (other_sq[dose].size()) {
                        Eligibility_Group* other_eg = other_sq[dose].top();
                        Eligibility_Group* eg = new Eligibility_Group();
                        eg->eligibility_day = other_eg->eligibility_day;
                        for (int bin : vac_campaign->get_unique_age_bins()) {
                            for (Person* p : other_eg->eligible_people[bin]) {
                                eg->eligible_people[bin].push_back(person_ptr_map[p]);
                            }
                        }
                        sq[dose].push(eg);
                        other_sq[dose].pop();
                    }

                    while (other_uq[dose].size()) {
                        Eligibility_Group* other_eg = other_uq[dose].top();
                        Eligibility_Group* eg = new Eligibility_Group();
                        eg->eligibility_day = other_eg->eligibility_day;
                        for (int bin : vac_campaign->get_unique_age_bins()) {
                            for (Person* p : other_eg->eligible_people[bin]) {
                                eg->eligible_people[bin].push_back(person_ptr_map[p]);
                            }
                        }
                        uq[dose].push(eg);
                        other_uq[dose].pop();
                    }
                }

                // TODO: copy doses available and pointers to them
                vac_campaign->copy_doses_available(o.vac_campaign);

                vac_campaign->set_potential_vaccinees(STANDARD_QUEUE, psv);
                vac_campaign->set_potential_vaccinees(URGENT_QUEUE, puv);
                vac_campaign->set_eligibility_queue(STANDARD_QUEUE, sq);
                vac_campaign->set_eligibility_queue(URGENT_QUEUE, uq);
            }
        }

        ~Community();

        Date* get_date();
        bool loadPopulation(std::string populationFilename, std::string comorbidityFilename = "", std::string publicActivityFilename = "", std::string immunityFilename = "");
        bool loadLocations(std::string locationFilename, std::string networkFilename = "");
        size_t getNumPeople() const { return _people.size(); }
        size_t getNumLocations() const { return _location.size(); }
        std::vector<Person*> getPeople() const { return _people; }
        size_t getNumInfected(int day); // includes people in incubation period
        size_t getNumInfectious(int day);
        size_t getNumSymptomatic(int day);
        size_t getNumNaive();
        Person* getPersonByID(int id);
        Location* getLocationByID(int id);
        Infection*  infect(int id, StrainType strain);
        int getDay() { return _day; }                                // what day is it?
        //void swapImmuneStates();
        void updatePersonStatus();
        // void tallyInfectionByLoc(Infection* inf);
        void updateHotLocations();
        void tick();                                                   // simulate one day

        void within_household_transmission();
        void between_household_transmission();
        void workplace_transmission();
        void school_transmission();
        void hospital_transmission();
        void nursinghome_transmission();
        void public_activity();
        void clear_public_activity();

        double social_distancing(int);
        void vaccinate();
        //void targetVaccination(Person* p); // routine vaccination on target birthday
        void updateVaccination();          // for boosting and multi-does vaccines
        void setVES(double f);
        void setVESs(std::vector<double> f);
        void setVac_Campaign(Vac_Campaign* vc) { vac_campaign = vc; }
        Vac_Campaign* getVac_Campaign() const { return vac_campaign; }
        std::vector<size_t> getNumNewlyInfected() { return cmty_ledger->_numNewlyInfected; }
        std::vector<size_t> getNumNewInfections(StrainType strain) { return cmty_ledger->_numNewInfectionsByStrain.at(strain); }
        std::vector<size_t> getNumNewInfectionsByLoc(string key) { return cmty_ledger->_numNewlyInfectedByLoc[key]; }
        std::vector<size_t> getNumNewlySymptomatic() { return cmty_ledger->_numNewlySymptomatic; }
        std::vector<size_t> getNumNewlySevere() { return cmty_ledger->_numNewlySevere; }
        std::vector<size_t> getNumNewlyCritical() { return cmty_ledger->_numNewlyCritical; }
        std::vector<size_t> getNumNewlyDead() { return cmty_ledger->_numNewlyDead; }

        std::vector<size_t> getNumVaccinatedCases() { return cmty_ledger->_numVaccinatedCases; }
        std::vector<size_t> getNumSeverePrev() { return cmty_ledger->_numSeverePrev; }
        std::vector<size_t> getNumHospInc() { return cmty_ledger->_numHospInc; }
        std::vector<size_t> getNumHospPrev() { return cmty_ledger->_numHospPrev; }
        std::vector<size_t> getNumIcuInc() { return cmty_ledger->_numIcuInc; }
        std::vector<size_t> getNumIcuPrev() { return cmty_ledger->_numIcuPrev; }
        std::vector<size_t> getNumDetectedCasesOnset() { return cmty_ledger->_numDetectedCasesOnset; }
        std::vector<size_t> getNumDetectedCasesReport() { return cmty_ledger->_numDetectedCasesReport; }
        std::vector<size_t> getNumDetectedHospitalizations() { return cmty_ledger->_numDetectedHospitalizations; }
        //std::vector<size_t> getNumDetectedDeaths() { return _numDetectedDeaths; }
        std::vector<size_t> getNumDetectedDeathsOnset() { return cmty_ledger->_numDetectedDeathsOnset; }
        std::vector<size_t> getNumDetectedDeathsReport() { return cmty_ledger->_numDetectedDeathsReport; }
        std::vector<pair<size_t, double>> getMeanNumSecondaryInfections() const ;
        std::vector<size_t> getCumulIncidenceByOutcome() { return cmty_ledger->_cumulIncByOutcome; }
        size_t getCumulIncidenceByOutcome( OutcomeType ot ) { return cmty_ledger->_cumulIncByOutcome[ot]; }

        double doSerosurvey (const ImmuneStateType ist, const std::vector<Person*> &pop, int time);
        double getHouseholdSecondaryAttackRate(std::vector<Person*> &pop);

        void flagInfectedLocation(Person* person, double relInfectiousness, LocationType locType, Location* _pLoc, int day);
        Infection* trace_contact(Person* &infecter, Location* source_loc, const map<double, vector<Person*>> &infectious_groups); //TODO: rename to findInfector()

        void reportCase(int onsetDate, long int reportDate, bool hospitalized);
        void reportDeath(int eventDate, long int reportDate);
        void tallyOutcome(OutcomeType ot) { cmty_ledger->_cumulIncByOutcome[ot]++; }
        vector< set<Person*, PerPtrComp> > traceForwardContacts();


//        int ageIntervalSize(int ageMin, int ageMax) { return std::accumulate(_personAgeCohortSizes+ageMin, _personAgeCohortSizes+ageMax,0); }

        void reset();                                                // reset the state of the community
        const std::vector<Location*> getLocations() const { return _location; }
        std::set<Location*, Location::LocPtrComp> getLocationsByType(LocationType lt) const { return _location_map.at(lt); }
        const std::vector<Person*> getAgeCohort(unsigned int age) const { assert(age<_personAgeCohort.size()); return _personAgeCohort[age]; }
        std::vector<double> getTimedIntervention(TimedIntervention ti) const { return cmty_ledger->_timedInterventions.at(ti); }
        void updateTimedIntervention(TimedIntervention ti, size_t date, double val) {
            const size_t current_size = cmty_ledger->_timedInterventions.at(ti).size();
            cmty_ledger->_timedInterventions.at(ti).resize(date);              // truncate
            cmty_ledger->_timedInterventions.at(ti).resize(current_size, val); // new values
        }
        void setSocialDistancingTimedIntervention(const vector<double> &vec) {
            cmty_ledger->_timedInterventions[SOCIAL_DISTANCING].clear();
            cmty_ledger->_timedInterventions[SOCIAL_DISTANCING] = vec;
            const double last_sd_value = cmty_ledger->_timedInterventions[SOCIAL_DISTANCING].back();
            cmty_ledger->_timedInterventions[SOCIAL_DISTANCING].resize(_par->runLength, last_sd_value);
        }
        void setSocialDistancingTimedIntervention(const vector<TimeSeriesAnchorPoint> &vec) {
            cmty_ledger->_timedInterventions[SOCIAL_DISTANCING].clear();
            cmty_ledger->_timedInterventions[SOCIAL_DISTANCING] = Date::linInterpolateTimeSeries(vec, _par->startJulianYear, _par->startDayOfYear);
            const double last_sd_value = cmty_ledger->_timedInterventions[SOCIAL_DISTANCING].back();
            cmty_ledger->_timedInterventions[SOCIAL_DISTANCING].resize(_par->runLength, last_sd_value);
        }
        void _clearSocialDistancingTimedIntervention() {
            cmty_ledger->_timedInterventions[SOCIAL_DISTANCING].clear();
        }
        void _extendSocialDistancingTimedIntervention(double val) {
            const size_t current_size = cmty_ledger->_timedInterventions.at(SOCIAL_DISTANCING).size();
            cmty_ledger->_timedInterventions[SOCIAL_DISTANCING].resize(current_size + _par->tuning_window, val);
        }
        void _reviseSocialDistancingTimedIntervention(double val) {
            const size_t current_size = cmty_ledger->_timedInterventions.at(SOCIAL_DISTANCING).size();
            assert(current_size >= _par->tuning_window);
            cmty_ledger->_timedInterventions[SOCIAL_DISTANCING].resize(current_size - _par->tuning_window);
            _extendSocialDistancingTimedIntervention(val);
        }
        double getTimedIntervention(TimedIntervention ti, size_t day) const { return cmty_ledger->_timedInterventions.at(ti)[day]; }

        vector<Location*> locsAtPixel(std::pair<double, double> px) { return _pixelMap[{px.first, px.second}]; }
        vector<Location*> locsAtPixel(double xP, double yP) { return _pixelMap[{xP, yP}]; }

        map<std::string, double> calculate_vax_stats(int day);
        vector<size_t> generateOffspringDistribution();

        CommunityLedger* get_ledger() const { return cmty_ledger; }
        void load_from_cache(CommunityLedger* cache_ledger, Date* cache_date, std::map<int, std::vector<Person*>> cache_hosp_ppl, Vac_Campaign* cache_vc);

        double getNumPeopleQuarantining(int day) {
            int tot = 0;
            for (Person* p : _people) { tot += (int) p->isQuarantining(day); }
            return (double) tot/_people.size();
        }

    protected:
        static const Parameters* _par;
        Date* _date;
        std::vector<Person*> _people;                                          // the array index is equal to the ID
        std::vector< std::vector<Person*> > _personAgeCohort;                  // array of pointers to people of the same age
        //int _personAgeCohortSizes[NUM_AGE_CLASSES];                          // size of each age cohort
        std::vector<Location*> _location;                                      // index is equal to the ID
        std::vector<Location*> _public_locations;                              // index is arbitrary
        std::map<LocationType, std::set<Location*, Location::LocPtrComp>> _location_map; //
        std::map<std::pair<double, double>, std::vector<Location*>> _pixelMap;
        //std::vector< std::vector<Person*> > _exposedQueue;                     // queue of people with n days of latency left
        int _day;                                                              // current day
        CommunityLedger* cmty_ledger;
        // std::vector<size_t> _numNewlyInfected;
        // std::map<StrainType, std::vector<size_t>> _numNewInfectionsByStrain;
        // std::map<std::string, std::vector<size_t>> _numNewlyInfectedByLoc;
        // std::vector<size_t> _numNewlySymptomatic;                              // true cases, no lag due to detection
        // std::vector<size_t> _numNewlySevere;                                   // true cases, no lag due to detection
        // std::vector<size_t> _numNewlyCritical;                                 // true cases, no lag due to detection
        // std::vector<size_t> _numNewlyDead;                                     // true cases, no lag due to detection
        // std::vector<size_t> _numVaccinatedCases;
        // std::vector<size_t> _numSeverePrev;
        // std::vector<size_t> _numHospInc;
        // std::vector<size_t> _numHospPrev;
        // std::vector<size_t> _numIcuInc;
        // std::vector<size_t> _numIcuPrev;
        // std::vector<size_t> _numDetectedCasesOnset;
        // std::vector<size_t> _numDetectedCasesReport;
        // std::vector<size_t> _numDetectedHospitalizations;
        // //std::vector<size_t> _numDetectedDeaths;
        // std::vector<size_t> _numDetectedDeathsOnset;
        // std::vector<size_t> _numDetectedDeathsReport;
        // std::vector<size_t> _cumulIncByOutcome;

        // groups of infectious people: indexed by day, location type, location ptr, and relative infectiousness
        // std::vector< std::map<LocationType, std::map<Location*, std::map<double, std::vector<Person*>>, Location::LocPtrComp>>> _isHot;
        // std::vector<Person*> _peopleByAge;
        Vac_Campaign* vac_campaign;
        // std::set<Person*> _revaccinate_set;          // not automatically re-vaccinated, just checked for boosting, multiple doses

        // std::map<TimedIntervention, std::vector<double>> timedInterventions;

        //bool _uniformSwap;                                            // use original swapping (==true); or parse swap file (==false)
        void _transmission(Location* source_loc, vector<Person*> at_risk_group, const map<double, vector<Person*>> &infectious_groups, const double T); // generic helper function

        void expandExposedQueues();
//        void _advanceTimers();
//        void _processBirthday(Person* p);
        //void _swapIfNeitherInfected(Person* p, Person* donor);
};
#endif
