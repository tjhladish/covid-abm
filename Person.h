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

struct Detection {
    OutcomeType detected_state;
    int reported_time;
};

class Infection {
    friend class Person;
    Infection() {
        infectedBegin     = INT_MAX;
        infectiousBegin   = INT_MAX;
        symptomBegin      = INT_MAX;
        severeBegin       = INT_MAX;
        hospitalizedBegin = INT_MAX;
        criticalBegin     = INT_MAX;
        icuBegin          = INT_MAX;
        deathTime         = INT_MAX;

        infectedPlace     = INT_MIN;
        infectedByID      = INT_MIN;
        infectiousEnd     = INT_MIN;
        symptomEnd        = INT_MIN;
        severeEnd         = INT_MIN;
        criticalEnd       = INT_MIN;

        relInfectiousness = 1.0;
        _detection        = nullptr;
    };

    ~Infection() {
        if (_detection) delete _detection;
    }

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
    std::vector<Infection*> infections_caused;
    double relInfectiousness;
    Detection* _detection;

    static const Parameters* _par;

  public:
    bool isLocallyAcquired()      const { return infectedByID != -1; }

    int getInfectedTime()         const { return infectedBegin; }
    int getInfectiousTime()       const { return infectiousBegin; }
    int getSymptomTime()          const { return symptomBegin; }
    int getSymptomEndTime()       const { return symptomEnd; }
    int getSevereTime()           const { return severeBegin; }
    int getHospitalizedTime()     const { return hospitalizedBegin; }
    int getCriticalTime()         const { return criticalBegin; }
    int getIcuTime()              const { return icuBegin; }
    int getDeathTime()            const { return deathTime; }
    double getRelInfectiousness() const { return relInfectiousness; }            // Could be time-varying in the future

    bool infected()               const { return infectedBegin     != INT_MAX; } // These functions check whether these things
    bool infectious()             const { return infectiousBegin   != INT_MAX; } // happen at any point during this infection, e.g.
    bool symptomatic()            const { return symptomBegin      != INT_MAX; } // have these been modified from their defaults?
    bool severe()                 const { return severeBegin       != INT_MAX; }
    bool critical()               const { return criticalBegin     != INT_MAX; }
    bool hospital()               const { return hospitalizedBegin != INT_MAX; }
    bool icu()                    const { return icuBegin          != INT_MAX; }
    bool fatal()                  const { return deathTime         != INT_MAX; }

    // if we ensure that death coincides with the end of symptoms/severity/criticality, then we don't also need to check deathtime
    bool isInfected(int now)      const { return infectedBegin <= now     and now < infectiousEnd;}
    bool isInfectious(int now)    const { return infectiousBegin <= now   and now < infectiousEnd;}
    bool isSymptomatic(int now)   const { return symptomBegin <= now      and now < symptomEnd; }
    bool isSevere(int now)        const { return severeBegin <= now       and now < severeEnd; }
    bool isCritical(int now)      const { return criticalBegin <= now     and now < criticalEnd; }
    bool inHospital(int now)      const { return hospitalizedBegin <= now and now < severeEnd;} // if fatal, severe state ends upon death
    bool inIcu(int now)           const { return icuBegin <= now          and now < criticalEnd;}
    bool isDead(int now)          const { return deathTime <= now; }

    void detect(OutcomeType detected_state, int report_date) { _detection = new Detection{detected_state, report_date}; }
    Detection* detection() { return _detection; }
    void log_transmission(Infection* inf) { infections_caused.push_back(inf); }
    size_t secondary_infection_tally () const { return infections_caused.size(); }
    std::vector<Infection*> get_infections_caused() { return infections_caused; }

    std::vector<int> generation_times () const {
        std::vector<int> times;
        for (Infection* other: infections_caused) { times.push_back(other->getInfectedTime() - getInfectedTime()); }
        return times;
    }

    void dumper() const {
        cerr << "Infection attributes:" << endl;
        cerr << "\tinfectedBegin    : " << infectedBegin            << endl;
        cerr << "\tinfectedPlace    : " << infectedPlace            << endl;
        cerr << "\tinfectedByID     : " << infectedByID             << endl;
        cerr << "\tinfectiousBegin  : " << infectiousBegin          << endl;
        cerr << "\tinfectiousEnd    : " << infectiousEnd            << endl;
        cerr << "\tsymptomBegin     : " << symptomBegin             << endl;
        cerr << "\tsymptomEnd       : " << symptomEnd               << endl;
        cerr << "\tsevereBegin      : " << severeBegin              << endl;
        cerr << "\tsevereEnd        : " << severeEnd                << endl;
        cerr << "\thospitalizedBegin: " << hospitalizedBegin        << endl;
        cerr << "\tcriticalBegin    : " << criticalBegin            << endl;
        cerr << "\tcriticalEnd      : " << criticalEnd              << endl;
        cerr << "\ticuBegin         : " << icuBegin                 << endl;
        cerr << "\tdeathTime        : " << deathTime                << endl;
        cerr << "\tinfections_caused: " << infections_caused.size() << endl;
        cerr << "\trelInfectiousness: " << relInfectiousness        << endl;
    }
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

        size_t getDaysImmune() const { return daysImmune; }
        void setDaysImmune(size_t di) { daysImmune = di; }

        bool getLongTermCare() { return long_term_care; }
        void setLongTermCare(bool b) { long_term_care = b; }

        ComorbidType hasComorbidity() { return comorbidity; }
        void setComorbidity(ComorbidType status) { comorbidity = status; }

        Location* getHomeLoc() { return home_loc; }
        void setHomeLoc(Location* loc) { home_loc = loc; }
        Location* getDayLoc() { return day_loc; }
        void setDayLoc(Location* loc) { day_loc = loc; }

        Location* getHospital() const { return home_loc->getHospital(); }
        void goToHospital() { getHospital()->addPerson(this); }
        void leaveHospital() { getHospital()->removePerson(this); }
//        void setImmunity() { immune = true; }
//        void copyImmunity(const Person *p);
        void resetImmunity();

        bool isCrossProtected(int time) const; // assumes binary cross-immunity
        bool isVaccineProtected(int time) const;

        inline int getInfectedByID(int infectionsago=0)      const { return getInfection(infectionsago)->infectedByID; }
        inline int getInfectedPlace(int infectionsago=0)     const { return getInfection(infectionsago)->infectedPlace; }
        inline int getInfectedTime(int infectionsago=0)      const { return getInfection(infectionsago)->infectedBegin; }

        inline int getInfectiousTime(int infectionsago=0)    const { return getInfection(infectionsago)->infectiousBegin; }
        inline int getSymptomTime(int infectionsago=0)       const { return getInfection(infectionsago)->symptomBegin; }
        inline int getSevereTime(int infectionsago=0)        const { return getInfection(infectionsago)->severeBegin; }
        inline int getCriticalTime(int infectionsago=0)      const { return getInfection(infectionsago)->criticalBegin; }

        inline int getHospitalizedTime(int infectionsago=0)  const { return getInfection(infectionsago)->hospitalizedBegin; }
        inline int getIcuTime(int infectionsago=0)           const { return getInfection(infectionsago)->icuBegin; }
        inline int getDeathTime(int infectionsago=0)         const { return getInfection(infectionsago)->deathTime; }

        Infection* getInfection(int infectionsago=0)         const { return infectionHistory[getNumNaturalInfections() - 1 - infectionsago]; }
        double     getRelInfectiousness(int infectionsago=0) const { return getInfection(infectionsago)->relInfectiousness; }
        inline int getNumNaturalInfections()                 const { return infectionHistory.size(); }

        int getNumVaccinations()                             const { return vaccineHistory.size(); }
        const std::vector<int>& getVaccinationHistory()      const { return vaccineHistory; }
        const std::vector<Infection*>& getInfectionHistory() const { return infectionHistory; }
        int daysSinceVaccination(int time)                   const { assert( vaccineHistory.size() > 0); return time - vaccineHistory.back(); } // isVaccinated() should be called first
        double vaccineProtection(const int time) const;

        Infection* infect(int sourceid, const Date* date, int sourceloc);
        void processDeath(Infection &infection, const int time);
        inline bool infect(const Date* date) {return infect(INT_MIN, date, INT_MIN);}

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
        bool isNewlyDead(int time)      const { return infectionHistory.size() > 0 and time == infectionHistory.back()->deathTime; }

        bool isVaccinated() const { return vaccineHistory.size() > 0; } // has been vaccinated
        bool isInfectable(int time) const;
        double remainingEfficacy(const int time) const;

        bool isNaive() const { return immune_state == NAIVE; }
                                                                        // does this person's immune state permit vaccination?
                                                                        // NB: inaccurate test results are possible
        bool isSeroEligible(VaccineSeroConstraint vsc, double falsePos, double falseNeg) const;
        bool vaccinate(int time);                                       // vaccinate this person
        static void setPar(const Parameters* par) { _par = par; }

        Infection& initializeNewInfection();
        Infection& initializeNewInfection(int time, size_t incubation_period, int sourceloc, int sourceid);

        static void reset_ID_counter() { NEXT_ID = 0; }
        bool isSurveilledPerson() { return id < _par->numSurveilledPeople; }

        void dumper() const {
            cerr << "Person ID: " << id << endl;
            cerr << "\thome loc: " << home_loc->getID() << endl;
            cerr << "\tday loc: " << (day_loc ? day_loc->getID() : -1) << endl;
            cerr << "\tday loc type: " << (day_loc ? day_loc->getType() : NUM_OF_LOCATION_TYPES) << endl;
            cerr << "\thospital: " << getHospital()->getID() << endl;
            cerr << "\tage: " << age << endl;
            cerr << "\tsex: " << sex << endl;
            cerr << "\tltc: " << long_term_care << endl;
            cerr << "\tinfection history size: " << infectionHistory.size() << endl;
        }


    protected:
        size_t id;                                                      // unique identifier
        Location* home_loc;                                             // family membership
        Location* day_loc;                                              // ID of location of work
        int age;                                                        // age in years
        SexType sex;                                                    // sex (gender)
        bool long_term_care;                                            // resident of nursing home, etc.
        ComorbidType comorbidity;                                       // person has a relevant comorbidity
        bool naiveVaccineProtection;                                    // if vaccinated, do we use the naive or non-naive VE_S?

        ImmuneStateType immune_state;
        std::vector<Infection*> infectionHistory;
        size_t daysImmune;                                              // number of days this person retains natural immunity
        std::vector<int> vaccineHistory;
        void clearInfectionHistory();

        static const Parameters* _par;
        static size_t NEXT_ID;                                          // unique ID to assign to the next Person allocated
};
#endif
