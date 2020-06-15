// Person.cpp

#include <cstdlib>
#include <cstring>
#include <climits>

#include <iostream>
#include <string>
#include <math.h>
#include <algorithm>

#include <assert.h>
#include <bitset>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Person.h"
#include "Community.h"
#include "Parameters.h"

using namespace covid::standard;

size_t Person::NEXT_ID = 0;

const Parameters* Person::_par;

Person::Person() {
    id = NEXT_ID++;
    age = -1;
    home_loc = nullptr;
    day_loc = nullptr;
//    home_id = -1;
//    day_id = -1;
//    for(int i=0; i<(int) NUM_OF_TIME_PERIODS; i++) _pLocation[i] = NULL;
    immune_state = NAIVE;
    //dead = false;
    //vaccinated = false;
    naiveVaccineProtection = false;
}


Person::~Person() {
    clearInfectionHistory();
}

void Person::clearInfectionHistory() {
    for (unsigned int i = 0; i < infectionHistory.size(); i++) {
        delete infectionHistory[i];
    }
    infectionHistory.clear();
}

Infection& Person::initializeNewInfection() {
//    setImmunity();
    Infection* infection = new Infection();
    infectionHistory.push_back(infection);

    switch( immune_state ) {
        case NAIVE:
        case NATURAL:
            immune_state = NATURAL;
            break;
        case VACCINATED:
        case NATURAL_AND_VACCINATED:
            immune_state = NATURAL_AND_VACCINATED;
            break;
        default:
            cerr << "ERROR: Unknown immune_state: " << immune_state << endl;
            exit(-1);
            break;
    }

    return *infection;
}


Infection& Person::initializeNewInfection(int time, int sourceloc, int sourceid) {
    Infection& infection      = initializeNewInfection();
    infection.infectedBegin   = time;
    infection.infectedPlace   = sourceloc;
    infection.infectedByID    = sourceid;
    infection.infectiousBegin = time + INFECTIOUSNESS_ONSET;
    infection.infectiousEnd   = infection.infectiousBegin + INFECTIOUS_PERIOD;
    return infection;
}

/*
// copyImmunity - copy immune status from person* p
void Person::copyImmunity(const Person* p) {
    assert(p!=NULL);
    immune = p->immune;
    _bVaccinated = p->_bVaccinated;

    vaccineHistory.clear();
    vaccineHistory.assign(p->vaccineHistory.begin(), p->vaccineHistory.end());

    clearInfectionHistory();
    for (int i=0; i < p->getNumNaturalInfections(); i++) {
        infectionHistory.push_back( new Infection(p->infectionHistory[i]) );
    }
}*/


// resetImmunity - reset immune status (infants)
void Person::resetImmunity() {
   // immune = false;
    clearInfectionHistory();
    vaccineHistory.clear();
    naiveVaccineProtection = false;
    immune_state = NAIVE;
}

/*
bool Person::naturalDeath(int t) {
    if (_nLifespan<=_nAge+(t/365.0)) {
        dead = true;
        return true;
    }
    return false;
}*/


bool Person::isInfectable(int time) const {
    return gsl_rng_uniform(RNG) < _par->susceptibilityByAge[age]
              and (isNaive() or (!isCrossProtected(time) and !isVaccineProtected(time)));
}


double Person::remainingEfficacy(const int time) const {
    double remainingFraction = 1.0;
    if (not isVaccinated()) {
        remainingFraction = 0.0;
    } else {
        if (_par->linearlyWaningVaccine) {
            // reduce by fraction of immunity duration that has waned
            int time_since_vac = daysSinceVaccination(time);
            if (time_since_vac > _par->vaccineImmunityDuration) {
                remainingFraction = 0.0;
            } else {
                remainingFraction -= ((double) time_since_vac) / _par->vaccineImmunityDuration;
            }
        }
    }
    return remainingFraction;
}


double Person::vaccineProtection(const int time) const {
    double ves;
    if (not isVaccinated()) {
        ves = 0.0;
    } else {
        if (daysSinceVaccination(time) > _par->vaccineImmunityDuration) {
            ves = 0.0;
        } else {
            if (naiveVaccineProtection == true) {
                ves = _par->VES_NAIVE;
            } else {
                ves = _par->VES;
            }
            ves *= remainingEfficacy(time);
        }
    }
    return ves;
}

//enum MaternalEffect { MATERNAL_PROTECTION, NO_EFFECT, MATERNAL_ENHANCEMENT };

/*MaternalEffect _maternal_antibody_effect(Person* p, const Parameters* _par, int time) {
    MaternalEffect effect = NO_EFFECT;
    if (p->getAge() == 0 and time >= 0) {                       // this is an infant, and we aren't reloading an infection history
        Person* mom = p->getLocation(HOME_NIGHT)->findMom();    // find a cohabitating female of reproductive age
        if (mom and mom->getImmunityBitset().any()) {           // if there is one and she has an infection history
            if (gsl_rng_uniform(RNG) < _par->infantImmuneProb) {
                effect = MATERNAL_PROTECTION;
            } else if (gsl_rng_uniform(RNG) < _par->infantSevereProb) {
                effect = MATERNAL_ENHANCEMENT;
            }
        }
    }

    return effect;
}*/


void Person::processDeath(Infection &infection, const int deathTime) {
    infection.deathTime         = deathTime;
    infection.infectiousEnd     = min(infection.infectiousEnd, deathTime);
    infection.symptomEnd        = min(infection.symptomEnd, deathTime);
    infection.severeEnd         = min(infection.severeEnd, deathTime);
    infection.hospitalizedBegin = min(infection.hospitalizedBegin, deathTime);
    infection.criticalEnd       = min(infection.criticalEnd, deathTime);
}



// infect - infect this individual
// returns true if infection occurs
bool Person::infect(int sourceid, int time, int sourceloc) {
    // Bail now if this person can not become infected
    // Not quite the same as "susceptible"--this person may be e.g. partially immune
    // due to natural infection or vaccination
    if (not isInfectable(time)) return false;
    const double remaining_efficacy = remainingEfficacy(time);  // due to vaccination; needs to be called before initializing new infection (still true?)

    // Create a new infection record
    Infection& infection = initializeNewInfection(time, sourceloc, sourceid);

    double symptomatic_probability = _par->pathogenicityByAge[age];             // may be modified by vaccination
    const double severe_given_case = SEVERE_FRACTION;           // might become age-, sex- or co-morbidity-structured in the future
    const double critical_given_severe = CRITICAL_FRACTION;

    const double effective_VEP = isVaccinated() ? _par->VEP*remaining_efficacy : 0.0;        // reduced symptoms due to vaccine
    symptomatic_probability *= (1.0 - effective_VEP);
    assert(symptomatic_probability >= 0.0);
    assert(symptomatic_probability <= 1.0);

    // determine disease outcome and timings
    if ( gsl_rng_uniform(RNG) < symptomatic_probability ) {
//cerr << "symptomatic ";
        // This is a case
        infection.symptomBegin = time + SYMPTOM_ONSET;
        if ( not (gsl_rng_uniform(RNG) < severe_given_case) ) {
//cerr << "not severe ";
            // It does not become severe
            infection.symptomEnd = infection.symptomBegin + SYMPTOM_DURATION_MILD;
        } else {
//cerr << "severe ";
            // It does progress and become severe
            infection.severeBegin = infection.symptomBegin + PRE_SEVERE_SYMPTOMATIC;

            // Is this person hospitalized when their severe symptoms begin?
            const float hosp_prob = long_term_care ? LTC_SEVERE_TO_HOSPITAL : SEVERE_TO_HOSPITAL;
            bool hosp = false;
            if (gsl_rng_uniform(RNG) < hosp_prob) {
//cerr << "hospitalized ";
                // This person is hospitalized when symptoms become severe
                hosp = true;
                infection.hospitalizedBegin = infection.severeBegin;
            }

            if (not (gsl_rng_uniform(RNG) < critical_given_severe)) {
//cerr << "not critical ";
                // Severe, but does not become critical
                infection.severeEnd     = infection.severeBegin   + SEVERE_DURATION;
                infection.symptomEnd    = infection.symptomBegin  + SYMPTOM_DURATION_SEVERE; // total symptomatic period for severe cases
            } else {
//cerr << "critical ";
                // It does progress to critical disease
                infection.criticalBegin = infection.severeBegin   + PRE_CRITICAL_SEVERE;
                infection.criticalEnd   = infection.criticalBegin + CRITICAL_DURATION;
                infection.severeEnd     = infection.severeBegin   + SEVERE_DURATION_CRITICAL;
                infection.symptomEnd    = infection.symptomBegin  + SYMPTOM_DURATION_CRITICAL;

                // Is this person put in intensive care (and thus also hospitalized, if not previously)?
                const float icu_prob = hosp ? CRITICAL_TO_ICU_IF_HOSP : CRITICAL_TO_ICU_IF_NOT_HOSP;
                bool death = false;
                if (gsl_rng_uniform(RNG) < icu_prob) {
//cerr << "icu ";
                    // Patient goes to intensive care
                    infection.icuBegin = infection.criticalBegin;
                    if (not hosp) { infection.hospitalizedBegin = infection.icuBegin; } // if they weren't hospitalized before, they are now
                    death = gsl_rng_uniform(RNG) < ICU_CRITICAL_MORTALITY;
                    if (death) {
//cerr << "death ";
                        // uniform randomly chose a day from the critical duration when death happens
                        processDeath(infection, infection.criticalBegin + gsl_rng_uniform_int(RNG, infection.criticalEnd - infection.criticalBegin));
                    }
                } else {
                    death = gsl_rng_uniform(RNG) < NON_ICU_CRITICAL_MORTALITY;
                    if (death) {
//cerr << "death ";
                        // death happens when critical symptoms begin, as this person is not receiving care
                        processDeath(infection, infection.criticalBegin);
                    }
                }
            }
/*            if (not isVaccinated() or gsl_rng_uniform(RNG) > _par->VEH*remaining_efficacy) { // Is this person unvaccinated or vaccinated but unlucky?
                infection.recoveryTime = infection.infectiousBegin + INFECTIOUS_PERIOD_SEVERE;
            }*/
        }
//cerr << endl;
    } else {
//cerr << "asymptomatic\n";
    }

    // Flag locations with (non-historical) infections, so that we know to look there for human->mosquito transmission
    // Negative days are historical (pre-simulation) events, and thus we don't care about modeling transmission
    for (int day = std::max(infection.infectiousBegin, 0); day < infection.infectiousEnd; day++) {
        Community::flagInfectedLocation(getHomeLoc(), day);
    }

    // if the antibody-primed vaccine-induced immunity can be acquired retroactively, upgrade this person from naive to mature
    if (_par->retroactiveMatureVaccine) naiveVaccineProtection = false;

    return true;
}


/*
bool Person::isWithdrawn(int time) const {
    if (infectionHistory.size() > 0) {
        Infection* infection = infectionHistory.back();
        if (time >= infection->withdrawnTime and time < infection->recoveryTime and not dead) {
            return true;
        }
    }
    return false;
}*/


bool Person::isCrossProtected(int time) const { // assumes binary cross-immunity
    return (getNumNaturalInfections() > 0) and  // has any past infection
           (infectionHistory.back()->infectedBegin + _par->daysImmune > time); // prev. infection w/in crossprotection period
}


bool Person::isVaccineProtected(int time) const {
    return isVaccinated() and
           ( !_par->vaccineLeaky or // if the vaccine isn't leaky
           (gsl_rng_uniform(RNG) < vaccineProtection(time)) ); // or it protects (i.e., doesn't leak this time)
}


bool Person::isSeroEligible(VaccineSeroConstraint vsc, double falsePos, double falseNeg) const {
    if (vsc == VACCINATE_ALL_SERO_STATUSES) return true;

    assert(falsePos >= 0.0 and falsePos <= 1.0);
    assert(falseNeg >= 0.0 and falseNeg <= 1.0);
    assert(vsc == VACCINATE_SEROPOSITIVE_ONLY or vsc == VACCINATE_SERONEGATIVE_ONLY);

    bool isSeroPos = not isNaive(); // fully susceptible == seronegative == false

    if ((isSeroPos and (falseNeg > gsl_rng_uniform(RNG)))       // sero+ but tests negative
        or (!isSeroPos and (falsePos > gsl_rng_uniform(RNG)))) { // sero- but tests positive
        isSeroPos = !isSeroPos;
    }

    bool eligible = vsc == VACCINATE_SEROPOSITIVE_ONLY ? isSeroPos : not isSeroPos;

    return eligible;
}


bool Person::vaccinate(int time) {
    if (!isDead(time)) {
        vaccineHistory.push_back(time);

        if ( isNaive() ) {
            naiveVaccineProtection = true;
        } else {
            naiveVaccineProtection = false;
        }

        if ( _par->vaccineLeaky == false ) { // all-or-none VE_S protection
            if ( (isNaive() and gsl_rng_uniform(RNG) < _par->VES_NAIVE) // vac someone who's naive
                 or gsl_rng_uniform(RNG) < _par->VES ) {                          // vac someone previously infected
                 switch( immune_state ) {
                    case NAIVE:
                    case VACCINATED:
                        immune_state = VACCINATED;
                        break;
                    case NATURAL:
                    case NATURAL_AND_VACCINATED:
                        immune_state = NATURAL_AND_VACCINATED;
                        break;
                    default:
                        cerr << "ERROR: Unknown immune_state: " << immune_state << endl;
                        exit(-1);
                        break;
                }

            }
        }
        return true;
    } else {
        return false;
    }
}
