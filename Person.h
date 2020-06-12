// Person.h
// A single individual.

#ifndef __PERSON_H
#define __PERSON_H
#include <bitset>
#include <vector>
#include <climits>
#include "Parameters.h"
#include "Location.h"

class Location;

class Infection {
    friend class Person;
    Infection() {
        infectedBegin     = INT_MIN;
        infectedPlace     = INT_MIN;
        infectedByID      = INT_MIN;
        infectiousBegin   = INT_MIN;
        infectiousEnd     = INT_MIN;

        symptomBegin      = INT_MIN;
        symptomEnd        = INT_MIN;

        severeBegin       = INT_MIN;
        severeEnd         = INT_MIN;
        hospitalizedBegin = INT_MIN;

        criticalBegin     = INT_MIN;
        criticalEnd       = INT_MIN;
        icuBegin          = INT_MIN;

        deathTime         = INT_MAX;        // this ONE should default to INT_MAX, the others INT_MIN
    };

    int infectedBegin;                      // when infected?
    int infectedPlace;                      // where infected?
    int infectedByID;                       // who infected this person
    int infectiousBegin;                    // when infectious period starts
    int infectiousEnd;

    int symptomBegin;                       // when symptoms start
    int symptomEnd;

    int severeBegin;
    int severeEnd;
    int hospitalizedBegin;

    int criticalBegin;
    int criticalEnd;
    int icuBegin;

    int deathTime;

  public:
    bool isLocallyAcquired()    const { return infectedByID != -1; }

    int getInfectedTime()       const { return infectedBegin; }
    int getInfectiousTime()     const { return infectiousBegin; }
    int getSymptomTime()        const { return symptomBegin; }
    int getSevereTime()         const { return severeBegin; }
    int getHospitalizedTime()   const { return hospitalizedBegin; }
    int getCriticalTime()       const { return criticalBegin; }
    int getIcuTime()            const { return icuBegin; }
    int getDeathTime()          const { return deathTime; }

//    int getInfectiousEndTime()  const { return infectiousDuration; }
    // if we ensure that death coincides with the end of symptoms/severity/criticality, then we don't also need to check deathtime
    bool isInfected(int now)    const { return infectedBegin <= now     and now < infectiousEnd;}
    bool isInfectious(int now)  const { return infectiousBegin <= now   and now < infectiousEnd;}
    bool isSymptomatic(int now) const { return symptomBegin <= now      and now < symptomEnd; }
    bool isSevere(int now)      const { return severeBegin <= now       and now < severeEnd; }
    bool isCritical(int now)    const { return criticalBegin <= now     and now < criticalEnd; }
    bool inHospital(int now)    const { return hospitalizedBegin <= now and now < severeEnd;}
    bool inIcu(int now)         const { return icuBegin <= now          and now < criticalEnd;}
    bool isDead(int now)        const { return deathTime <= now; }
};

class Person {
    public:
        Person();
        ~Person();
        inline int getID() const { return id; }

        int getAge() const { return age; }
        void setAge(int a) { age = a; }

        SexType getSex() const { return sex; }
        void setSex(SexType s) { sex = s; }

        bool getLongTermCare() { return long_term_care; }
        void setLongTermCare(bool b) { long_term_care = b; }

        Location* getHomeLoc() { return home_loc; }
        void setHomeLoc(Location* loc) { home_loc = loc; }
        Location* getDayLoc() { return day_loc; }
        void setDayLoc(Location* loc) { day_loc = loc; }

//        int getHomeID() const { return home_id; }
//        void setHomeID(int n) { home_id = n; }
//        int getDayID() const { return day_id; }
//        void setDayID(int n) { day_id = n; }
//        void setImmunity() { immune = true; }
//        void copyImmunity(const Person *p);
        void resetImmunity();

//        bool isSusceptible() const;                  // is susceptible (and is alive)
        bool isCrossProtected(int time) const; // assumes binary cross-immunity
        bool isVaccineProtected(int time) const;

//        inline Location* getLocation(TimePeriod timeofday) const { return _pLocation[(int) timeofday]; }
//        inline void setLocation(Location* p, TimePeriod timeofday) { _pLocation[(int) timeofday] = p; }

        inline int getInfectedByID(int infectionsago=0)     const { return getInfection(infectionsago)->infectedByID; }
        inline int getInfectedPlace(int infectionsago=0)    const { return getInfection(infectionsago)->infectedPlace; }
        inline int getInfectedTime(int infectionsago=0)     const { return getInfection(infectionsago)->infectedBegin; }

        inline int getInfectiousTime(int infectionsago=0)   const { return getInfection(infectionsago)->infectiousBegin; }
//        inline int getInfectiousDuration(int infectionsago=0)   const { return getInfection(infectionsago)->infectiousDuration; }

        inline int getSymptomTime(int infectionsago=0)      const { return getInfection(infectionsago)->symptomBegin; }
//        inline int getSymptomDuration(int infectionsago=0)      const { return getInfection(infectionsago)->symptomDuration; }

        inline int getSevereTime(int infectionsago=0)       const { return getInfection(infectionsago)->severeBegin; }
//        inline int getSevereDuration(int infectionsago=0)   const { return getInfection(infectionsago)->severeDuration; }

        inline int getCriticalTime(int infectionsago=0)     const { return getInfection(infectionsago)->criticalBegin; }
//        inline int getCriticalDuration(int infectionsago=0) const { return getInfection(infectionsago)->criticalDuration; }

        inline int getHospitalizationTime(int infectionsago=0) const { return getInfection(infectionsago)->hospitalizedBegin; }
        inline int getIcuTime(int infectionsago=0) const { return getInfection(infectionsago)->icuBegin; }
        inline int getDeathTime(int infectionsago=0)        const { return getInfection(infectionsago)->deathTime; }

        const Infection* getInfection(int infectionsago=0) const { return infectionHistory[getNumNaturalInfections() - 1 - infectionsago]; }
        //inline void setRecoveryTime(int time, int infectionsago=0) { infectionHistory[getNumNaturalInfections() - 1 - infectionsago]->recoveryTime = time; }
        inline int getNumNaturalInfections() const { return infectionHistory.size(); }

        int getNumVaccinations() const { return vaccineHistory.size(); }
        const std::vector<int>& getVaccinationHistory() const { return vaccineHistory; }
        const std::vector<Infection*>& getInfectionHistory() const { return infectionHistory; }
        int daysSinceVaccination(int time) const { assert( vaccineHistory.size() > 0); return time - vaccineHistory.back(); } // isVaccinated() should be called first
        double vaccineProtection(const int time) const;

        bool infect(int sourceid, int time, int sourceloc);
        void processDeath(Infection &infection, const int time);
        inline bool infect(int time) {return infect(INT_MIN, time, INT_MIN);}

        // TODO -- the following functions assume that only the most recent infection needs to be inspected
        bool inHospital(int time) const { return infectionHistory.size() > 0 and infectionHistory.back()->inHospital(time); }
        bool inIcu(int time)      const { return infectionHistory.size() > 0 and infectionHistory.back()->inIcu(time); }
        bool isDead(int time)     const { return infectionHistory.size() > 0 and infectionHistory.back()->isDead(time); } // no deaths due to other causes

        bool isNewlyInfected(int time)  const { return infectionHistory.size() > 0 and time == infectionHistory.back()->infectedBegin; }
        bool isInfected(int time)       const { return infectionHistory.size() > 0 and infectionHistory.back()->isInfected(time); }
        bool isInfectious(int time)     const { return infectionHistory.size() > 0 and infectionHistory.back()->isInfectious(time); }
        bool isSymptomatic(int time)    const { return infectionHistory.size() > 0 and infectionHistory.back()->isSymptomatic(time); }
        bool isSevere(int time)         const { return infectionHistory.size() > 0 and infectionHistory.back()->isSevere(time); }
        bool isCritical(int time)       const { return infectionHistory.size() > 0 and infectionHistory.back()->isCritical(time); }



//        bool isNewlyInfected(int time) const;                           // became infected today?
//        bool isInfectious(int time) const;
//        bool isInfected(int time) const;                                // either exposed or infectious
//        bool isSymptomatic(int time) const;                             // has symptoms
//        bool isSevere(int time) const;                                  // used for estimating hospitalizations
//        bool isCritical(int time) const;                                // used for estimating ICU demand
        bool isVaccinated() const { return vaccineHistory.size() > 0; } // has been vaccinated
        bool isInfectable(int time) const;
        double remainingEfficacy(const int time) const;

        bool isNaive() const { return immune_state == NAIVE; }
                                                                      // does this person's immune state permit vaccination?
                                                                      // NB: inaccurate test results are possible
        bool isSeroEligible(VaccineSeroConstraint vsc, double falsePos, double falseNeg) const;
        bool vaccinate(int time);                                     // vaccinate this person
        static void setPar(const Parameters* par) { _par = par; }

        //static const double _fIncubationDistribution[MAX_INCUBATION];

        Infection& initializeNewInfection();
        Infection& initializeNewInfection(int time, int sourceloc, int sourceid);

        static void reset_ID_counter() { NEXT_ID = 1; }

    protected:
        size_t id;                                                  // unique identifier
        Location* home_loc;                                                // family membership
        Location* day_loc;                                                 // ID of location of work
//        Location *_pLocation[(int) NUM_OF_TIME_PERIODS];            // where this person is at morning, day, and evening
        int age;                                                    // age in years
        SexType sex;                                                // sex (gender)
        bool long_term_care;                                        // resident of nursing home, etc.
//        bool dead;                                                  // is dead
//        bool vaccinated;                                            // has been vaccinated
        bool naiveVaccineProtection; // if vaccinated, do we use the naive or non-naive VE_S?

        ImmuneStateType immune_state;
        std::vector<Infection*> infectionHistory;
        std::vector<int> vaccineHistory;
        void clearInfectionHistory();

        static const Parameters* _par;
        static size_t NEXT_ID;                                          // unique ID to assign to the next Person allocated
};
#endif
