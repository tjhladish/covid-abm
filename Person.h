// Person.h
// A single individual.

#ifndef __PERSON_H
#define __PERSON_H
#include <bitset>
#include <vector>
#include <climits>
#include "Parameters.h"
#include "Location.h"

enum QuarantineLevel {FULL, MODERATE, MINIMAL, NUM_OF_QUARANTINE_LEVELS};

class Location;
class Community;
class Vaccinee;

struct Detection {
    Detection() {
        detected_state = NUM_OF_OUTCOME_TYPES;
        reported_time = INT_MAX;
    };

    Detection(OutcomeType ds, int rt) : detected_state(ds), reported_time(rt) {};

    Detection(const Detection& o) {
        detected_state = o.detected_state;
        reported_time  = o.reported_time;
    };

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

        infectedPlace     = nullptr;
        infectedBy        = nullptr;
        infectiousEnd     = INT_MIN;
        symptomEnd        = INT_MIN;
        severeEnd         = INT_MIN;
        criticalEnd       = INT_MIN;

        strain            = WILDTYPE;
        relInfectiousness = 1.0;
        _detection        = nullptr;
    };

    Infection(const Infection& o) {
        infectedBegin     = o.infectedBegin;
        infectedPlace     = o.infectedPlace;
        infectedBy        = o.infectedBy;
        infectiousBegin   = o.infectiousBegin;
        infectiousEnd     = o.infectiousEnd;
        symptomBegin      = o.symptomBegin;
        symptomEnd        = o.symptomEnd;
        severeBegin       = o.severeBegin;
        severeEnd         = o.severeEnd;
        hospitalizedBegin = o.hospitalizedBegin;
        criticalBegin     = o.criticalBegin;
        criticalEnd       = o.criticalEnd;
        icuBegin          = o.icuBegin;
        deathTime         = o.deathTime;
        infections_caused = o.infections_caused;
        strain            = o.strain;
        relInfectiousness = o.relInfectiousness;
        _detection        = o._detection ? new Detection(*(o._detection)) : nullptr;
    };

    ~Infection() {
        if (_detection) delete _detection;
    }

    int infectedBegin;                      // when infected?
    Location* infectedPlace;                // where infected?
    Person* infectedBy;                     // who infected this person
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
    StrainType strain;
    double relInfectiousness;
    Detection* _detection;

    static const Parameters* _par;

  public:
    bool isLocallyAcquired()      const { return infectedBy != nullptr; }

    int getInfectedTime()         const { return infectedBegin; }
    int getInfectiousTime()       const { return infectiousBegin; }
    int getSymptomTime()          const { return symptomBegin; }
    int getSymptomEndTime()       const { return symptomEnd; }
    int getSevereTime()           const { return severeBegin; }
    int getHospitalizedTime()     const { return hospitalizedBegin; }
    int getCriticalTime()         const { return criticalBegin; }
    int getIcuTime()              const { return icuBegin; }
    int getDeathTime()            const { return deathTime; }
    StrainType getStrain()        const { return strain; }
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

    Detection* getDetection()        { return _detection; }
    bool isDetected()          const { return (bool) _detection; }
    bool isDetectedOn(int now) const { return (bool) isDetected() and _detection->reported_time == now; }

    void log_transmission(Infection* inf) { infections_caused.push_back(inf); }
    size_t secondary_infection_tally () const { return infections_caused.size(); }
    std::vector<Infection*> get_infections_caused() { return infections_caused; }

    std::vector<int> generation_times () const {
        std::vector<int> times;
        for (Infection* other: infections_caused) { times.push_back(other->getInfectedTime() - getInfectedTime()); }
        return times;
    }

    void dumper() const;
};

class Person {
    friend Vaccinee;
    public:
        Person();
        Person(const Person& o) {
            id                      = o.id;
            home_loc                = o.home_loc;
            day_loc                 = o.day_loc;
            patronized_locs         = o.patronized_locs;
            age                     = o.age;
            sex                     = o.sex;
            long_term_care          = o.long_term_care;
            comorbidity             = o.comorbidity;
            naiveVaccineProtection  = o.naiveVaccineProtection;
            immune_state            = o.immune_state;
            infectionHistory        = std::vector<Infection*>(o.infectionHistory.size());
            quarantineStart         = o.quarantineStart;
            quarantineEnd           = o.quarantineEnd;

            for(size_t i = 0; i < o.infectionHistory.size(); ++i) {
                infectionHistory[i] = new Infection(*(o.infectionHistory[i]));
            }

            naturalImmunityDuration = o.naturalImmunityDuration;
            startingNaturalEfficacy = o.startingNaturalEfficacy;
            immunityQuantile        = o.immunityQuantile;
            vaccineHistory          = o.vaccineHistory;
        };

        ~Person();

        struct PerPtrComp { bool operator()(const Person* A, const Person* B) const { return A->getID() < B->getID(); } };

        inline int getID() const { return id; }

        int getAge() const { return age; }
        void setAge(int a) { age = a; }

        SexType getSex() const { return sex; }
        void setSex(SexType s) { sex = s; }

        bool isHCW() { return getDayLoc() and ((getDayLoc()->getType() == HOSPITAL) or (getDayLoc()->getType() == NURSINGHOME)); }

        double getStartingNaturalEfficacy() const { return startingNaturalEfficacy; }
        void setStartingNaturalEfficacy(double sne) { startingNaturalEfficacy = sne; }
        //size_t getDaysImmune() const { return daysImmune; }
        //void setDaysImmune(size_t di) { daysImmune = di; }

        int getNaturalImmunityDuration() const { return naturalImmunityDuration; }
        void setNaturalImmunityDuration(int nid) { naturalImmunityDuration = nid; }

        int getVaccineImmunityDuration(size_t dose, StrainType strain) const { return _par->immunityDuration(immunityQuantile, dose, strain); }

        double getImmunityQuantile() const { return immunityQuantile; }
        void setImmunityQuantile(double iq) { immunityQuantile = iq; }

        bool getLongTermCare() { return long_term_care; }
        void setLongTermCare(bool b) { long_term_care = b; }

        ComorbidType hasComorbidity() { return comorbidity; }
        void setComorbidity(ComorbidType status) { comorbidity = status; }

        Location* getHomeLoc() { return home_loc; }
        void setHomeLoc(Location* loc) { home_loc = loc; }
        double getRiskiness() const { return home_loc->getRiskiness(); }
        Location* getDayLoc() { return day_loc; }
        void setDayLoc(Location* loc) { day_loc = loc; }

        Location* getHospital() const { return home_loc->getHospital(); }
        void goToHospital() { getHospital()->addPerson(this); }
        void leaveHospital() { getHospital()->removePerson(this); }

        void addPatronizedLocation(Location* loc) { patronized_locs.push_back(loc); } // person is sometimes a customer of these businesses
        vector<Location*> getPatronizedLocations() const { return patronized_locs; }
        void setPatronizedLocations(vector<Location*> pl) { patronized_locs = pl; }
//        void setImmunity() { immune = true; }
//        void copyImmunity(const Person *p);
        void resetImmunity();

        bool isCrossProtected(int time, StrainType strain) const; // assumes binary cross-immunity
        bool isVaccineProtected(int time, StrainType strain) const;

        inline Person* getInfectedBy(int infectionsago=0)      const { return getInfection(infectionsago)->infectedBy; }
        inline Location* getInfectedPlace(int infectionsago=0) const { return getInfection(infectionsago)->infectedPlace; }
        inline int getInfectedTime(int infectionsago=0)        const { return getInfection(infectionsago)->infectedBegin; }

        inline int getInfectiousTime(int infectionsago=0)      const { return getInfection(infectionsago)->infectiousBegin; }
        inline int getSymptomTime(int infectionsago=0)         const { return getInfection(infectionsago)->symptomBegin; }
        inline int getSevereTime(int infectionsago=0)          const { return getInfection(infectionsago)->severeBegin; }
        inline int getCriticalTime(int infectionsago=0)        const { return getInfection(infectionsago)->criticalBegin; }

        inline int getHospitalizedTime(int infectionsago=0)    const { return getInfection(infectionsago)->hospitalizedBegin; }
        inline int getIcuTime(int infectionsago=0)             const { return getInfection(infectionsago)->icuBegin; }
        inline int getDeathTime(int infectionsago=0)           const { return getInfection(infectionsago)->deathTime; }

        Infection* getInfection(int infectionsago=0)           const { return infectionHistory[getNumNaturalInfections() - 1 - infectionsago]; }
        StrainType getStrain(int infectionsago=0)              const { return getInfection(infectionsago)->strain; }
        double     getRelInfectiousness(int infectionsago=0)   const { return getInfection(infectionsago)->relInfectiousness; }
        inline int getNumNaturalInfections()                   const { return infectionHistory.size(); }

        int getNumVaccinations()                               const { return vaccineHistory.size(); }
        const std::vector<int>& getVaccinationHistory()        const { return vaccineHistory; }
        const std::vector<Infection*>& getInfectionHistory()   const { return infectionHistory; }
        int daysSinceVaccination(int time)                     const { assert( vaccineHistory.size() > 0); return time - vaccineHistory.back(); } // isVaccinated() should be called first
        double vaccineProtection(const int time, const StrainType strain) const;

        // strain determined by source, unless source is nullptr
        Infection* infect(Community* community, Person* source, const Date* date, Location* sourceloc, StrainType strain = NUM_OF_STRAIN_TYPES, bool check_susceptibility = true);
        void processDeath(Community* community, Infection &infection, const int time);
        inline Infection* infect(Community* community, const Date* date, StrainType strain) {return infect(community, nullptr, date, nullptr, strain);}

        // TODO -- the following functions assume that only the most recent infection needs to be inspected
        bool inHospital(int time) const { return infectionHistory.size() > 0 and infectionHistory.back()->inHospital(time); }
        bool inIcu(int time)      const { return infectionHistory.size() > 0 and infectionHistory.back()->inIcu(time); }
        bool isDead(int time)     const { return infectionHistory.size() > 0 and infectionHistory.back()->isDead(time); } // no deaths due to other causes
        bool isAlive(int time)    const { return not isDead(time); }

        bool isNewlyInfected(int time)  const { return infectionHistory.size() > 0 and time == infectionHistory.back()->infectedBegin; }
        bool isInfected(int time)       const { return infectionHistory.size() > 0 and infectionHistory.back()->isInfected(time); }
        bool isInfectious(int time)     const { return infectionHistory.size() > 0 and infectionHistory.back()->isInfectious(time); }
        bool isSymptomatic(int time)    const { return infectionHistory.size() > 0 and infectionHistory.back()->isSymptomatic(time); }
        bool isSevere(int time)         const { return infectionHistory.size() > 0 and infectionHistory.back()->isSevere(time); }
        bool isCritical(int time)       const { return infectionHistory.size() > 0 and infectionHistory.back()->isCritical(time); }
        bool isNewlyDead(int time)      const { return infectionHistory.size() > 0 and time == infectionHistory.back()->deathTime; }

        bool isVaccinated() const { return vaccineHistory.size() > 0; } // has been vaccinated
        bool isInfectable(int time, StrainType strain) const;
        double remainingEfficacy(const int time) const;

        bool hasBeenInfected() const { return (bool) infectionHistory.size(); }
        //bool isNaive() const { return immune_state == NAIVE; }
        bool isImmuneState(ImmuneStateType is) const { return is == immune_state; }

                                                                        // does this person's immune state permit vaccination?
                                                                        // NB: inaccurate test results are possible
        //bool isSeroEligible(VaccineSeroConstraint vsc, double falsePos, double falseNeg) const;
        bool isSeroEligible() const;
        static void setPar(const Parameters* par) { _par = par; }

        Infection& initializeNewInfection();
        Infection& initializeNewInfection(int time, size_t incubation_period, Location* sourceloc, Person* source);

        static void reset_ID_counter() { NEXT_ID = 0; }
        bool isSurveilledPerson() { return id < _par->numSurveilledPeople; }

        void selfQuarantine(const size_t today, const size_t quarantineDuration);
        void endQuarantine();
        bool isQuarantining(const size_t today);

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
            cerr << "\tinfection days: "; for (size_t i = 0; i < infectionHistory.size(); ++i) { cerr << getInfectedTime(i) << ' '; } cerr << endl;
            cerr << "\tvaccination history size: " << vaccineHistory.size() << endl;
            cerr << "\tvaccine dose days: "; for (int d : vaccineHistory) { cerr << d << ' '; } cerr << endl;
        }


    protected:
        size_t id;                                                      // unique identifier
        Location* home_loc;                                             // family membership
        Location* day_loc;                                              // ID of location of work
        std::vector<Location*> patronized_locs;                         // business patronized by this person
        int age;                                                        // age in years
        SexType sex;                                                    // sex (gender)
        bool long_term_care;                                            // resident of nursing home, etc.
        ComorbidType comorbidity;                                       // person has a relevant comorbidity
        bool naiveVaccineProtection;                                    // if vaccinated, do we use the naive or non-naive VE_S?

        ImmuneStateType immune_state;
        std::vector<Infection*> infectionHistory;

        int naturalImmunityDuration;                                    // number of days this person retains natural immunity (for all-or-none immunity model)
        double startingNaturalEfficacy;                                 // level of protection (Pr{reisting infection}, analogous to VE_S) acquired by person after natural infection
        double immunityQuantile;                                        // represents the individual's immune durability

        std::vector<int> vaccineHistory;                                // vector of days on which vaccinations were received
        void clearInfectionHistory();

        size_t quarantineStart;
        size_t quarantineEnd;

        bool vaccinate(int time);                                       // vaccinate this person

        static const Parameters* _par;
        static size_t NEXT_ID;                                          // unique ID to assign to the next Person allocated
};
#endif
