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

class Community {
    public:
        Community(const Parameters* parameters, Date* date);
        Community(const Community& o) {
            _cumulIncByOutcome           = o._cumulIncByOutcome;// is static
            _par                         = o._par;// is static const
            _date                        = o._date;                  // is static
            _people                      = std::vector<Person*>(o._people.size());
            _personAgeCohort             = o._personAgeCohort;
            _location                    = std::vector<Location*>(o._location.size());
            _public_locations            = o._public_locations;
            _location_map                = o._location_map;
            _exposedQueue                = o._exposedQueue;
            _day                         = o._day;
            _numNewlyInfected            = o._numNewlyInfected;
            _numNewInfectionsByStrain    = o._numNewInfectionsByStrain;
            _numNewlySymptomatic         = o._numNewlySymptomatic;
            _numNewlySevere              = o._numNewlySevere;
            _numNewlyCritical            = o._numNewlyCritical;
            _numNewlyDead                = o._numNewlyDead;
            _numVaccinatedCases          = o._numVaccinatedCases;
            _numSeverePrev               = o._numSeverePrev;
            _numHospInc                  = o._numHospInc;
            _numHospPrev                 = o._numHospPrev;
            _numIcuInc                   = o._numIcuInc;
            _numIcuPrev                  = o._numIcuPrev;
            _numDetectedCasesOnset       = o._numDetectedCasesOnset;
            _numDetectedCasesReport      = o._numDetectedCasesReport;
            _numDetectedHospitalizations = o._numDetectedHospitalizations;
            _numDetectedDeaths           = o._numDetectedDeaths;
            _isHot                       = o._isHot;
            _peopleByAge                 = o._peopleByAge;
            vac_campaign                 = o.vac_campaign ? new Vac_Campaign(*(o.vac_campaign)) : nullptr;
            _revaccinate_set             = o._revaccinate_set;// is static
            timedInterventions           = o.timedInterventions;// is static

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
            }

            deque<Person*> urgent_queue = o.vac_campaign->getUrgentQueue();
            for(Person* &p : urgent_queue) { p = person_ptr_map[p]; }
            vac_campaign->setUrgentQueue(urgent_queue);

            deque<Person*> standard_queue = o.vac_campaign->getStandardQueue();
            for(Person* &p : standard_queue) { p = person_ptr_map[p]; }
            vac_campaign->setStandardQueue(standard_queue);

            vector< set<Person*> > revaccinate_queue;
            for(auto &q : o.vac_campaign->getRevaccinateQueue()) {
                set<Person*> daily_queue;
                for(Person* p : q) { daily_queue.insert(person_ptr_map[p]); }
                revaccinate_queue.push_back(daily_queue);
            }
            vac_campaign->setRevaccinateQueue(revaccinate_queue);

        }

        ~Community();

        Date* get_date();
        bool loadPopulation(std::string populationFilename, std::string comorbidityFilename = "", std::string publicActivityFilename = "", std::string immunityFilename = "");
        bool loadLocations(std::string locationFilename, std::string networkFilename = "");
        size_t getNumPeople() const { return _people.size(); }
        std::vector<Person*> getPeople() const { return _people; }
        size_t getNumInfected(int day); // includes people in incubation period
        size_t getNumInfectious(int day);
        size_t getNumSymptomatic(int day);
        size_t getNumNaive();
        Person* getPersonByID(int id);
        Infection*  infect(int id, StrainType strain);
        int getDay() { return _day; }                                // what day is it?
        //void swapImmuneStates();
        void updatePersonStatus();
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
        std::vector<size_t> getNumNewlyInfected() { return _numNewlyInfected; }
        std::vector<size_t> getNumNewInfections(StrainType strain) { return _numNewInfectionsByStrain.at(strain); }
        std::vector<size_t> getNumNewlySymptomatic() { return _numNewlySymptomatic; }
        std::vector<size_t> getNumNewlySevere() { return _numNewlySevere; }
        std::vector<size_t> getNumNewlyCritical() { return _numNewlyCritical; }
        std::vector<size_t> getNumNewlyDead() { return _numNewlyDead; }

        std::vector<size_t> getNumVaccinatedCases() { return _numVaccinatedCases; }
        std::vector<size_t> getNumSeverePrev() { return _numSeverePrev; }
        std::vector<size_t> getNumHospInc() { return _numHospInc; }
        std::vector<size_t> getNumHospPrev() { return _numHospPrev; }
        std::vector<size_t> getNumIcuInc() { return _numIcuInc; }
        std::vector<size_t> getNumIcuPrev() { return _numIcuPrev; }
        std::vector<size_t> getNumDetectedCasesOnset() { return _numDetectedCasesOnset; }
        std::vector<size_t> getNumDetectedCasesReport() { return _numDetectedCasesReport; }
        std::vector<size_t> getNumDetectedHospitalizations() { return _numDetectedHospitalizations; }
        std::vector<size_t> getNumDetectedDeaths() { return _numDetectedDeaths; }
        std::vector<pair<size_t, double>> getMeanNumSecondaryInfections() const ;

        void flagInfectedLocation(Person* person, double relInfectiousness, LocationType locType, Location* _pLoc, int day);
        Infection* trace_contact(Person* &infecter, Location* source_loc, const map<double, vector<Person*>> &infectious_groups);
        void reportCase(int onsetDate, long int reportDate, bool hospitalized);
        void reportDeath(int eventDate, long int reportDate);

//        int ageIntervalSize(int ageMin, int ageMax) { return std::accumulate(_personAgeCohortSizes+ageMin, _personAgeCohortSizes+ageMax,0); }

        void reset();                                                // reset the state of the community
        const std::vector<Location*> getLocations() const { return _location; }
        const std::vector<Person*> getAgeCohort(unsigned int age) const { assert(age<_personAgeCohort.size()); return _personAgeCohort[age]; }
        std::vector<double> getTimedIntervention(TimedIntervention ti) const { return timedInterventions.at(ti); }
        void updateTimedIntervention(TimedIntervention ti, size_t date, double val) {
            const size_t current_size = timedInterventions.at(ti).size();
            timedInterventions.at(ti).resize(date);              // truncate
            timedInterventions.at(ti).resize(current_size, val); // new values
        }
        double getTimedIntervention(TimedIntervention ti, size_t day) const { return timedInterventions.at(ti)[day]; }

        static std::vector<size_t> _cumulIncByOutcome;
    protected:
        static const Parameters* _par;
        static Date* _date;
        std::vector<Person*> _people;                                          // the array index is equal to the ID
        std::vector< std::vector<Person*> > _personAgeCohort;                  // array of pointers to people of the same age
        //int _personAgeCohortSizes[NUM_AGE_CLASSES];                          // size of each age cohort
        std::vector<Location*> _location;                                      // index is equal to the ID
        std::vector<Location*> _public_locations;                              // index is arbitrary
        std::map<LocationType, std::set<Location*, Location::LocPtrComp>> _location_map; //
        std::vector< std::vector<Person*> > _exposedQueue;                     // queue of people with n days of latency left
        int _day;                                                              // current day
        std::vector<size_t> _numNewlyInfected;
        std::map<StrainType, std::vector<size_t>> _numNewInfectionsByStrain;
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
        std::vector<size_t> _numDetectedDeaths;

        // groups of infectious people: indexed by day, location type, location ptr, and relative infectiousness
        std::vector< std::map<LocationType, std::map<Location*, std::map<double, std::vector<Person*>>, Location::LocPtrComp>>> _isHot;
        static std::vector<Person*> _peopleByAge;
        Vac_Campaign* vac_campaign;
        static std::set<Person*> _revaccinate_set;          // not automatically re-vaccinated, just checked for boosting, multiple doses
        std::map<TimedIntervention, std::vector<double>> timedInterventions;

        //bool _uniformSwap;                                            // use original swapping (==true); or parse swap file (==false)
        void _transmission(Location* source_loc, vector<Person*> at_risk_group, const map<double, vector<Person*>> &infectious_groups, const double T); // generic helper function

        void expandExposedQueues();
//        void _advanceTimers();
//        void _processBirthday(Person* p);
        //void _swapIfNeitherInfected(Person* p, Person* donor);
};
#endif
