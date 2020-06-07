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
        infectedTime     = INT_MIN;
        infectedPlace    = INT_MIN;
        infectedByID     = INT_MIN;
        infectiousTime   = INT_MIN;
        infectiousDuration = INT_MIN;

        symptomTime      = INT_MIN;
        symptomDuration  = INT_MIN;

        severeTime       = INT_MIN;
        severeDuration   = INT_MIN;
        hospitalizedTime = INT_MIN;

        criticalTime     = INT_MIN;
        criticalDuration = INT_MIN;
        icuTime          = INT_MIN;

        deathTime        = INT_MAX;
    };

    int infectedTime;                               // when infected?
    int infectedPlace;                              // where infected?
    int infectedByID;                               // who infected this person
    int infectiousTime;                             // when infectious period starts
    int infectiousDuration;                         //

    int symptomTime;                                // when symptoms start
    int symptomDuration;                            //

    int severeTime;
    int severeDuration;
    int hospitalizedTime;

    int criticalTime;
    int criticalDuration;
    int icuTime;

    int deathTime;

  public:
    // TODO - consider storing state end-points rather than durations.  likely fewer operations
    bool isLocallyAcquired()    const { return infectedByID != -1; }
    int getInfectedTime()       const { return infectedTime; }
    int getInfectiousTime()     const { return infectiousTime; }
    int getInfectiousDuration() const { return infectiousDuration; }
    // if we ensure that death coincides with the end of symptoms/severity/criticality, then we don't also need to check deathtime
    bool isInfected(int now)    const { return infectedTime <= now and now < (infectiousTime + infectiousDuration);}
    bool isInfectious(int now)  const { return infectiousTime <= now and now < (infectiousTime + infectiousDuration);}
    bool isSymptomatic(int now) const { return symptomTime <= now and now < (symptomTime + symptomDuration) and now < deathTime; }
    bool isSevere(int now)      const { return severeTime <= now and now < (severeTime + severeDuration); }
    bool isCritical(int now)    const { return criticalTime <= now and now < (criticalTime + criticalDuration); }
    bool isDead(int now)        const { return deathTime <= now; }
    bool inHospital(int now)    const { return hospitalizedTime <= now and now < (severeTime + severeDuration);}
    bool inICU(int now)         const { return icuTime <= now and now < (criticalTime + criticalDuration);}
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
        inline int getInfectedTime(int infectionsago=0)     const { return getInfection(infectionsago)->infectedTime; }

        inline int getInfectiousTime(int infectionsago=0)   const { return getInfection(infectionsago)->infectiousTime; }
        inline int getInfectiousDuration(int infectionsago=0)   const { return getInfection(infectionsago)->infectiousDuration; }

        inline int getSymptomTime(int infectionsago=0)      const { return getInfection(infectionsago)->symptomTime; }
        inline int getSymptomDuration(int infectionsago=0)      const { return getInfection(infectionsago)->symptomDuration; }

        inline int getSevereTime(int infectionsago=0)       const { return getInfection(infectionsago)->severeTime; }
        inline int getSevereDuration(int infectionsago=0)   const { return getInfection(infectionsago)->severeDuration; }

        inline int getCriticalTime(int infectionsago=0)     const { return getInfection(infectionsago)->criticalTime; }
        inline int getCriticalDuration(int infectionsago=0) const { return getInfection(infectionsago)->criticalDuration; }

        inline int getHospitalizationTime(int infectionsago=0) const { return getInfection(infectionsago)->hospitalizedTime; }
        inline int getIcuTime(int infectionsago=0) const { return getInfection(infectionsago)->icuTime; }
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
        inline bool infect(int time) {return infect(INT_MIN, time, INT_MIN);}

        bool inHospital(int time) const { return infectionHistory.size() > 0 and infectionHistory.back()->inHospital(time); }
        bool inICU(int time)      const { return infectionHistory.size() > 0 and infectionHistory.back()->inICU(time); }
        bool isDead(int time)     const { return infectionHistory.size() > 0 and infectionHistory.back()->isDead(time); } // no deaths due to other causes

        bool isNewlyInfected(int time)  const { return infectionHistory.size() > 0 and time == infectionHistory.back()->infectedTime; }
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
