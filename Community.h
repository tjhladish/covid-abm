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
        Community(const Parameters* parameters);
        virtual ~Community();
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
        bool infect(int id, int day);
        int getDay() { return _day; }                                // what day is it?
        //void swapImmuneStates();
        void updateDiseaseStatus();
        void tick(int day);                                           // simulate one day
        void within_household_transmission();
        void between_household_transmission();
        void workplace_and_school_transmission();
        void location_transmission(std::set<Location*, Location::LocPtrComp> &locations);

        //void setNoSecondaryTransmission() { _bNoSecondaryTransmission = true; }

        float social_distancing(int) const;
        void vaccinate(CatchupVaccinationEvent cve);
        void targetVaccination(Person* p); // routine vaccination on target birthday
        void updateVaccination();          // for boosting and multi-does vaccines
        void setVES(double f);
        void setVESs(std::vector<double> f);
        std::vector<size_t> getNumNewlyInfected() { return _numNewlyInfected; }
        std::vector<size_t> getNumNewlySymptomatic() { return _numNewlySymptomatic; }
        std::vector<size_t> getNumVaccinatedCases() { return _numVaccinatedCases; }
        std::vector<size_t> getNumSevereCases() { return _numSevereCases; }
        static void flagInfectedLocation(Location* _pLoc, int day);

//        int ageIntervalSize(int ageMin, int ageMax) { return std::accumulate(_personAgeCohortSizes+ageMin, _personAgeCohortSizes+ageMax,0); }

        void reset();                                                // reset the state of the community
        const std::vector<Location*> getLocations() const { return _location; }
        const std::vector<Person*> getAgeCohort(unsigned int age) const { assert(age<_personAgeCohort.size()); return _personAgeCohort[age]; }

    protected:
        static const Parameters* _par;
        std::vector<Person*> _people;                                          // the array index is equal to the ID
        std::vector< std::vector<Person*> > _personAgeCohort;                  // array of pointers to people of the same age
        //int _personAgeCohortSizes[NUM_AGE_CLASSES];                            // size of each age cohort
        //double *_fMortality;                                                 // mortality by year, starting from 0
        std::vector<Location*> _location;                                      // the array index is equal to the ID
        std::map<LocationType, std::set<Location*, Location::LocPtrComp>> _location_map; //
        std::vector< std::vector<Person*> > _exposedQueue;                     // queue of people with n days of latency left
        int _day;                                                              // current day
        std::vector<size_t> _numNewlyInfected;
        std::vector<size_t> _numNewlySymptomatic;
        std::vector<size_t> _numVaccinatedCases;
        std::vector<size_t> _numSevereCases;
        static std::vector<std::set<Location*, Location::LocPtrComp> > _isHot;
        static std::vector<Person*> _peopleByAge;
        static std::map<int, std::set<std::pair<Person*, Person*> > > _delayedBirthdays;
        static std::set<Person*> _revaccinate_set;          // not automatically re-vaccinated, just checked for boosting, multiple doses

        //bool _uniformSwap;                                            // use original swapping (==true); or parse swap file (==false)

        void expandExposedQueues();
//        void _advanceTimers();
//        void _processBirthday(Person* p);
//        void _processDelayedBirthdays();
        //void _swapIfNeitherInfected(Person* p, Person* donor);
};
#endif
