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
    home_id = -1;
    day_id = -1;
    for(int i=0; i<(int) NUM_OF_TIME_PERIODS; i++) _pLocation[i] = NULL;
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
    Infection& infection = initializeNewInfection();
    infection.infectedTime  = time;
    infection.infectedPlace = sourceloc;
    infection.infectedByID  = sourceid; // TODO - What kind of ID is this?
    //infection.infectiousTime = Parameters::sampler(INCUBATION_CDF, gsl_rng_uniform(RNG)) + time;
    infection.infectiousTime = time + _par->infectiousOnset;
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
    return isNaive() or (!isCrossProtected(time) and !isVaccineProtected(time));
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


// infect - infect this individual
// returns true if infection occurs
bool Person::infect(int sourceid, int time, int sourceloc) {
    // Bail now if this person can not become infected
    // Not quite the same as "susceptible"--this person may be e.g. partially immune
    // due to natural infection or vaccination
    if (not isInfectable(time)) return false;

    //const int numPrevInfections = getNumNaturalInfections();  // these both need to be called
    const double remaining_efficacy = remainingEfficacy(time);  // before initializing new infection

    // Create a new infection record
    Infection& infection = initializeNewInfection(time, sourceloc, sourceid);

    double symptomatic_probability = _par->basePathogenicity * SYMPTOMATIC_BY_AGE[age];
    const double severe_given_case = _par->severeFraction;
    // TODO -- add support for *all* outcomes

    if (symptomatic_probability > 1.0) symptomatic_probability = 1.0;
    const double effective_VEP = isVaccinated() ? _par->VEP*remaining_efficacy : 0.0;        // reduced symptoms due to vaccine
    symptomatic_probability *= (1.0 - effective_VEP);
    assert(symptomatic_probability >= 0.0);
    assert(symptomatic_probability <= 1.0);

    infection.recoveryTime = infection.infectiousTime + INFECTIOUS_PERIOD_ASYMPTOMATIC;            // may be changed below 

    if ( gsl_rng_uniform(RNG) < symptomatic_probability ) {         // Is this a case?
        const double severe_rand = gsl_rng_uniform(RNG);
        infection.recoveryTime = infection.infectiousTime + INFECTIOUS_PERIOD_MILD;                // may yet be changed below 
        if ( severe_rand < severe_given_case ) {                     // Is this a severe case?
            if (not isVaccinated() or gsl_rng_uniform(RNG) > _par->VEH*remaining_efficacy) { // Is this person unvaccinated or vaccinated but unlucky?
                infection.recoveryTime = infection.infectiousTime + INFECTIOUS_PERIOD_SEVERE;
//                infection.severeDisease = true; // TODO -- SET severeTime
            }
        }

        // Determine if this person withdraws (stops going to work/school)
/*        infection.symptomTime = infection.infectiousTime + SYMPTOMATIC_DELAY;
        const int symptomatic_duration = infection.recoveryTime - infection.symptomTime;
        const int symptomatic_active_period = gsl_ran_geometric(RNG, 0.5) - 1; // min generator value is 1 trial
        infection.withdrawnTime = symptomatic_active_period < symptomatic_duration ?
                                  infection.symptomTime + symptomatic_active_period :
                                  infection.withdrawnTime;*/
    }

    // Flag locations with (non-historical) infections, so that we know to look there for human->mosquito transmission
    // Negative days are historical (pre-simulation) events, and thus we don't care about modeling transmission
    for (int day = std::max(infection.infectiousTime, 0); day < infection.recoveryTime; day++) {
        for (int t=0; t<(int) NUM_OF_TIME_PERIODS; t++) {
            Community::flagInfectedLocation(_pLocation[t], day);
        }
    }

    // if the antibody-primed vaccine-induced immunity can be acquired retroactively, upgrade this person from naive to mature
    if (_par->retroactiveMatureVaccine) naiveVaccineProtection = false;

    return true;
}


bool Person::isNewlyInfected(int time) const {
    if (infectionHistory.size() > 0) {
        Infection* infection = infectionHistory.back();
        if (time == infection->infectedTime) {
            return true;
        }
    }
    return false;
}


bool Person::isInfected(int time) const {
    if (infectionHistory.size() > 0) {
        Infection* infection = infectionHistory.back();
        if (time >= infection->infectedTime and time < infection->recoveryTime) {
            return true;
        }
    }
    return false;
}


bool Person::isInfectious(int time) const { // TODO -- we will almost certainly want to make this more sophisticated
    if (infectionHistory.size() > 0) {
        Infection* infection = infectionHistory.back();
        if (time >= infection->getInfectiousTime() and time < (infection->getInfectiousTime() + _par->infectiousDuration) and not isDead(time)) {
            return true;
        }
    }
    return false;
}


bool Person::isSymptomatic(int time) const {
    if (infectionHistory.size() > 0) {
        return infectionHistory.back()->isSymptomatic(time);
    } else {
        return false;
    }
}


bool Person::isSevere(int time) const {
    if (infectionHistory.size() > 0) {
        return infectionHistory.back()->isSevere(time);
    } else {
        return false;
    }
}


bool Person::isCritical(int time) const {
    if (infectionHistory.size() > 0) {
        return infectionHistory.back()->isCritical(time);
    } else {
        return false;
    }
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
           (infectionHistory.back()->infectedTime + _par->daysImmune > time); // prev. infection w/in crossprotection period 
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
