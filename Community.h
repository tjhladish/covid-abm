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

class Person;

// We use this to created a vector of people, sorted by decreasing age.  Used for aging/immunity swapping.
//struct PerPtrComp { bool operator()(const Person* A, const Person* B) const { return A->getAge() > B->getAge(); } };

class Community {
    public:
        Community(const Parameters* parameters, Date* date);
        ~Community();

        Date* get_date();
        bool loadPopulation(std::string populationFilename, std::string immunityFilename = "");
        bool loadLocations(std::string locationFilename, std::string networkFilename = "");
        size_t getNumPeople() const { return _people.size(); }
        std::vector<Person*> getPeople() const { return _people; }
        size_t getNumInfected(int day); // includes people in incubation period
        size_t getNumInfectious(int day);
        size_t getNumSymptomatic(int day);
        size_t getNumNaive();
        //void populate(Person **parray, int targetpop);
        Person* getPersonByID(int id);
        bool infect(int id);
        int getDay() { return _day; }                                // what day is it?
        //void swapImmuneStates();
        void updateDiseaseStatus();
        void updateHotLocations();
        void tick();                                                   // simulate one day

        void within_household_transmission();
        void between_household_transmission();
        void workplace_transmission();
        void school_transmission();

        double social_distancing(int);
        void vaccinate(CatchupVaccinationEvent cve);
        void targetVaccination(Person* p); // routine vaccination on target birthday
        void updateVaccination();          // for boosting and multi-does vaccines
        void setVES(double f);
        void setVESs(std::vector<double> f);
        std::vector<size_t> getNumNewlyInfected() { return _numNewlyInfected; }
        std::vector<size_t> getNumNewlySymptomatic() { return _numNewlySymptomatic; }
        std::vector<size_t> getNumNewlyDead() { return _numNewlyDead; }
        std::vector<size_t> getNumVaccinatedCases() { return _numVaccinatedCases; }
        std::vector<size_t> getNumSeverePrev() { return _numSeverePrev; }
        std::vector<size_t> getNumHospInc() { return _numHospInc; }
        std::vector<size_t> getNumHospPrev() { return _numHospPrev; }
        std::vector<size_t> getNumIcuInc() { return _numIcuInc; }
        std::vector<size_t> getNumIcuPrev() { return _numIcuPrev; }
        std::vector<size_t> getNumDetectedCasesOnset() { return _numDetectedCasesOnset; }
        std::vector<size_t> getNumDetectedCasesReport() { return _numDetectedCasesReport; }
        std::vector<size_t> getNumDetectedDeaths() { return _numDetectedDeaths; }
        std::vector<pair<size_t, double>> getMeanNumSecondaryInfections() const ;

        static void flagInfectedLocation(LocationType locType, Location* _pLoc, int day);
        Infection* trace_contact(int &infectee_id, vector<Person*> &source_candidates, int infectious_count);
        static void reportCase(int onsetDate, long int reportDate);
        static void reportDeath(int eventDate, long int reportDate);

//        int ageIntervalSize(int ageMin, int ageMax) { return std::accumulate(_personAgeCohortSizes+ageMin, _personAgeCohortSizes+ageMax,0); }

        void reset();                                                // reset the state of the community
        const std::vector<Location*> getLocations() const { return _location; }
        const std::vector<Person*> getAgeCohort(unsigned int age) const { assert(age<_personAgeCohort.size()); return _personAgeCohort[age]; }
        std::vector<float> getTimedIntervention(TimedIntervention ti) const { return timedInterventions.at(ti); }
        void updateTimedIntervention(TimedIntervention ti, size_t date, float val) {
            const size_t current_size = timedInterventions.at(ti).size();
            timedInterventions.at(ti).resize(date);              // truncate
            timedInterventions.at(ti).resize(current_size, val); // new values
        }
        float getTimedIntervention(TimedIntervention ti, size_t day) const { return timedInterventions.at(ti)[day]; }

    protected:
        static const Parameters* _par;
        static Date* _date;
        std::vector<Person*> _people;                                          // the array index is equal to the ID
        std::vector< std::vector<Person*> > _personAgeCohort;                  // array of pointers to people of the same age
        //int _personAgeCohortSizes[NUM_AGE_CLASSES];                            // size of each age cohort
        //double *_fMortality;                                                 // mortality by year, starting from 0
        std::vector<Location*> _location;                                      // the array index is equal to the ID
        std::map<LocationType, std::set<Location*, Location::LocPtrComp>> _location_map; //
        std::vector< std::vector<Person*> > _exposedQueue;                     // queue of people with n days of latency left
        int _day;                                                              // current day
        std::vector<size_t> _numNewlyInfected;
        std::vector<size_t> _numNewlySymptomatic;                                // true cases, no lag
        std::vector<size_t> _numNewlyDead;                                       // true cases, no lag
        std::vector<size_t> _numVaccinatedCases;
        std::vector<size_t> _numSeverePrev;
        std::vector<size_t> _numHospInc;
        std::vector<size_t> _numHospPrev;
        std::vector<size_t> _numIcuInc;
        std::vector<size_t> _numIcuPrev;
        static std::vector<size_t> _numDetectedCasesOnset;
        static std::vector<size_t> _numDetectedCasesReport;
        static std::vector<size_t> _numDetectedDeaths;
        static std::vector<std::map<LocationType, std::map<Location*, int, Location::LocPtrComp>>> _isHot;
        static std::vector<Person*> _peopleByAge;
        static std::map<int, std::set<std::pair<Person*, Person*> > > _delayedBirthdays;
        static std::set<Person*> _revaccinate_set;          // not automatically re-vaccinated, just checked for boosting, multiple doses
        std::map<TimedIntervention, std::vector<float>> timedInterventions;

        //bool _uniformSwap;                                            // use original swapping (==true); or parse swap file (==false)
        void _transmission(Location* source_loc, vector<Person*> at_risk_group, const double T, const int infectious_count); // generic helper function

        void expandExposedQueues();
//        void _advanceTimers();
//        void _processBirthday(Person* p);
//        void _processDelayedBirthdays();
        //void _swapIfNeitherInfected(Person* p, Person* donor);
};
#endif
