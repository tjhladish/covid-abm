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
const Parameters* Infection::_par;

void Infection::dumper() const {
    cerr << "Infection attributes:" << endl;
    cerr << "\tstrain           : " << strain                                        << endl;
    cerr << "\tinfectedBegin    : " << infectedBegin                                 << endl;
    cerr << "\tinfectedPlaceID  : " << (infectedPlace ? infectedPlace->getID() : -1) << endl;
    cerr << "\tinfectedByID     : " << (infectedBy ? infectedBy->getID() : -1)       << endl;
    cerr << "\tinfectiousBegin  : " << infectiousBegin                               << endl;
    cerr << "\tinfectiousEnd    : " << infectiousEnd                                 << endl;
    cerr << "\tsymptomBegin     : " << symptomBegin                                  << endl;
    cerr << "\tsymptomEnd       : " << symptomEnd                                    << endl;
    cerr << "\tsevereBegin      : " << severeBegin                                   << endl;
    cerr << "\tsevereEnd        : " << severeEnd                                     << endl;
    cerr << "\thospitalizedBegin: " << hospitalizedBegin                             << endl;
    cerr << "\tcriticalBegin    : " << criticalBegin                                 << endl;
    cerr << "\tcriticalEnd      : " << criticalEnd                                   << endl;
    cerr << "\ticuBegin         : " << icuBegin                                      << endl;
    cerr << "\tdeathTime        : " << deathTime                                     << endl;
    cerr << "\tinfections_caused: " << infections_caused.size()                      << endl;
    cerr << "\trelInfectiousness: " << relInfectiousness                             << endl;
}


Person::Person() {
    id = NEXT_ID++;
    age = -1;
    home_loc = nullptr;
    day_loc = nullptr;
    immune_state = NAIVE;
    daysImmune = INT_MAX;
    naiveVaccineProtection = false;
    long_term_care = false;
    comorbidity = HEALTHY;
    quarantineStart = INT_MAX;
    quarantineEnd = INT_MAX;
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
    Community::_cumulIncByOutcome[ASYMPTOMATIC]++;

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


Infection& Person::initializeNewInfection(int time, size_t incubation_period, Location* sourceloc, Person* source) {
    Infection& infection      = initializeNewInfection();
    infection.infectedBegin   = time;
    infection.infectedPlace   = sourceloc;
    infection.infectedBy      = source;
    infection.strain          = source ? source->getStrain() : WILDTYPE;
    infection.infectiousBegin = time + _par->infectiousness_onset(incubation_period);
    infection.infectiousEnd   = time + incubation_period + SYMPTOMATIC_INFECTIOUS_PERIOD; // person may not actually be symptomatic!
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


bool Person::isInfectable(const int time, const StrainType strain) const {
    return (not isInfected(time)) and gsl_rng_uniform(RNG) < _par->susceptibilityByAge[age]
              and (isNaive() or (!isCrossProtected(time, strain) and !isVaccineProtected(time, strain)));
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


double Person::vaccineProtection(const int time, const StrainType strain) const {
    double ves;
    if (not isVaccinated()) {
        ves = 0.0;
    } else {
        const size_t dose = vaccineHistory.size() - 1;
        if (daysSinceVaccination(time) > _par->vaccineImmunityDuration) {
            ves = 0.0;
        } else {
            if (naiveVaccineProtection == true) {
                ves = _par->VES_NAIVE_at(dose, strain);
            } else {
                ves = _par->VES_at(dose, strain);
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
    Community::_cumulIncByOutcome[DEATH]++;
    infection.deathTime         = deathTime;
    infection.infectiousEnd     = min(infection.infectiousEnd, deathTime);
    infection.symptomEnd        = min(infection.symptomEnd, deathTime);
    infection.severeEnd         = min(infection.severeEnd, deathTime);
    infection.criticalEnd       = min(infection.criticalEnd, deathTime);
}



// infect - infect this individual
// returns non-null pointer if infection occurs
Infection* Person::infect(Person* source, const Date* date, Location* sourceloc, StrainType strain, bool check_susceptibility) {
    const int time = date->day();
    // Bail now if this person can not become infected
    // Not quite the same as "susceptible"--this person may be e.g. partially immune
    // due to natural infection or vaccination
    if (check_susceptibility and not isInfectable(time, strain)) { return nullptr; }
    const double remaining_efficacy = remainingEfficacy(time);  // due to vaccination; needs to be called before initializing new infection (still true?)

    // Create a new infection record
    const size_t incubation_period = _par->symptom_onset(); // may not be symptomatic, but this is used to determine infectiousness onset
    Infection& infection = initializeNewInfection(time, incubation_period, sourceloc, source);
    if (not source) { infection.strain = strain; }

    const size_t dose = vaccineHistory.size() - 1;
    const double effective_VEP = isVaccinated() ? _par->VEP_at(dose, strain)*remaining_efficacy : 0.0;        // reduced pathogenicity due to vaccine
    const double effective_VEH = isVaccinated() ? _par->VEH_at(dose, strain)*remaining_efficacy : 0.0;        // reduced severity (hospitalization) due to vaccine
    const double effective_VEF = isVaccinated() ? _par->VEF_at(dose, strain)*remaining_efficacy : 0.0;        // reduced fatality due to vaccine
    const double effective_VEI = isVaccinated() ? _par->VEI_at(dose, strain)*remaining_efficacy : 0.0;        // reduced infectiousness

    double symptomatic_probability = _par->pathogenicityByAge[age] * (1.0 - effective_VEP);           // may be modified by vaccination
    const double severe_given_case = _par->probSeriousOutcome.at(SEVERE)[comorbidity][age] * (1.0 - effective_VEH);
    const double critical_given_severe = _par->probSeriousOutcome.at(CRITICAL)[comorbidity][age];
    infection.relInfectiousness *= (1.0 - effective_VEI); // 0.25 b/c most (== asymptomatic) people are not very infectious
    if (infection.strain == B_1_1_7) {
        infection.relInfectiousness *= 1.6;
        symptomatic_probability *= 1.1;                  // TODO -- these values shouldn't be hard-coded here
    } else if (infection.strain == B_1_617_2) {
        infection.relInfectiousness *= 3.0;
        symptomatic_probability *= 1.2;
    }

    const double highly_infectious_threshold = 8.04; // 80th %ile for overall SARS-CoV-2 from doi: 10.7554/eLife.65774, "Fig 4-Fig Sup 3"

    // determine disease outcome and timings
    if ( gsl_rng_uniform(RNG) > symptomatic_probability ) {
       // asymptomatic; est mean relInfectiousness = 6.03
        infection.relInfectiousness *= gsl_ran_weibull(RNG, 6.72, 3.33) > highly_infectious_threshold ? 4.0 : 0.25;
    } else {
        // symptomatic; est mean relInfectiousness = 6.69
        infection.relInfectiousness *= gsl_ran_weibull(RNG, 7.40, 3.81) > highly_infectious_threshold ? 4.0 : 0.25;
        //const size_t symptom_onset = _par->symptom_onset();
        Community::_cumulIncByOutcome[MILD]++;
        infection.symptomBegin = time + incubation_period;
        if ( not (gsl_rng_uniform(RNG) < severe_given_case) ) {
            // It does not become severe
            infection.symptomEnd = infection.symptomBegin + _par->symptom_duration_mild();
        } else {
            // It does progress and become severe
            Community::_cumulIncByOutcome[SEVERE]++;
            infection.severeBegin = infection.symptomBegin + _par->pre_severe_symptomatic();

            // Is this person hospitalized when their severe symptoms begin?
            const float hosp_prob = long_term_care ? LTC_SEVERE_TO_HOSPITAL : SEVERE_TO_HOSPITAL;
            bool hosp = false;
            if (gsl_rng_uniform(RNG) < hosp_prob) {
                // This person is hospitalized when symptoms become severe
                hosp = true;
                infection.hospitalizedBegin = infection.severeBegin;
            }

            if (not (gsl_rng_uniform(RNG) < critical_given_severe)) {
                // Severe, but does not become critical
                infection.severeEnd     = infection.severeBegin   + _par->severe_only_duration();
                infection.symptomEnd    = infection.severeEnd; // TODO - extend symptoms beyond severe period (relevant for e.g. econ analyses)
            } else {
                // It does progress to critical disease
                Community::_cumulIncByOutcome[CRITICAL]++;
                const size_t severe_only_duration = _par->severe_only_duration();
                const size_t pre_critical_severe  = round((float) severe_only_duration / 2);
                const size_t post_critical_severe = severe_only_duration - pre_critical_severe;
                infection.criticalBegin = infection.severeBegin   + pre_critical_severe;
                infection.criticalEnd   = infection.criticalBegin + _par->critical_duration();
                infection.severeEnd     = infection.criticalEnd   + post_critical_severe;
                infection.symptomEnd    = infection.severeEnd; // TODO - extend symptoms beyond severe period (relevant for e.g. econ analyses)

                // Is this person put in intensive care (and thus also hospitalized, if not previously)?
                const float icu_prob = hosp ? CRITICAL_TO_ICU_IF_HOSP : CRITICAL_TO_ICU_IF_NOT_HOSP;
                bool death = false;
                if (gsl_rng_uniform(RNG) < icu_prob) {
                    // Patient goes to intensive care
                    infection.icuBegin = infection.criticalBegin;
                    if (not hosp) { infection.hospitalizedBegin = infection.icuBegin; } // if they weren't hospitalized before, they are now
                    death = gsl_rng_uniform(RNG) < _par->icuMortality(comorbidity, age, infection.icuBegin) * (1.0 - effective_VEF);
                    if (death) {
                        // uniform randomly chose a day from the critical duration when death happens
                        processDeath(infection, infection.criticalBegin + _par->sampleIcuTimeToDeath());
                    }
                } else {
                    death = gsl_rng_uniform(RNG) < NON_ICU_CRITICAL_MORTALITY * (1.0 - effective_VEF);
                    if (death) {
                        // non-icu death, happens when critical symptoms begin, as this person is not receiving care
                        processDeath(infection, infection.criticalBegin + _par->sampleCommunityTimeToDeath());
                    }
                }
            }
/*            if (not isVaccinated() or gsl_rng_uniform(RNG) > _par->VEH*remaining_efficacy) { // Is this person unvaccinated or vaccinated but unlucky?
                infection.recoveryTime = infection.infectiousBegin + INFECTIOUS_PERIOD_SEVERE;
            }*/
        }
    }

    // Detection/reporting!  TODO -- implement probabilistic detection in hospitals
    if (isSurveilledPerson() and time >= 0) {
        bool detected = false;
        OutcomeType detected_state = NUM_OF_OUTCOME_TYPES;
        int sample_collection_date = 0;
        long int report_date = 0;
        const size_t reporting_lag = _par->reportingLag(REPORTING_RNG, date);
        if (infection.infected() and gsl_rng_uniform(RNG) < _par->probFirstDetection[time][ASYMPTOMATIC]) {
            // extra delay, e.g. time during infection someone would be identified by chance screening
            const size_t infectious_period = infection.infectiousEnd - infection.infectiousBegin;
            const int tracing_lag = gsl_rng_uniform_int(RNG, infectious_period);
            sample_collection_date = infection.infectiousBegin + tracing_lag;
            report_date = sample_collection_date + reporting_lag;
            detected_state = ASYMPTOMATIC;
            detected = true;
        } else if (infection.symptomatic() and gsl_rng_uniform(RNG) < _par->probFirstDetection[time][MILD]) {
            sample_collection_date = infection.symptomBegin + _par->symptomToTestLag;
            report_date = sample_collection_date + reporting_lag;
            detected_state = MILD;
            detected = true;
        } else if (infection.severe() and (gsl_rng_uniform(RNG) < _par->probFirstDetection[time][SEVERE])) {
        //} else if (infection.severe() and (infection.inHospital(infection.severeBegin) or gsl_rng_uniform(RNG) < _par->probFirstDetection[time][SEVERE])) {
            sample_collection_date = infection.severeBegin;
            report_date = sample_collection_date + reporting_lag;
            detected_state = SEVERE;
            detected = true;
        } else if (infection.critical() and (gsl_rng_uniform(RNG) < _par->probFirstDetection[time][CRITICAL])) {
        //} else if (infection.critical() and (infection.inHospital(infection.criticalBegin) or gsl_rng_uniform(RNG) < _par->probFirstDetection[time][CRITICAL])) {
            sample_collection_date = infection.criticalBegin;
            report_date = sample_collection_date + reporting_lag;
            detected_state = CRITICAL;
            detected = true;
        } else if (infection.fatal() and gsl_rng_uniform(RNG) < _par->probFirstDetection[time][DEATH]) {
            sample_collection_date = infection.deathTime;
            report_date = sample_collection_date + _par->deathReportingLag(REPORTING_RNG);
            detected_state = DEATH;
            detected = true;
        }

        if (detected) {
            infection.detect(detected_state, report_date);
            Community::reportCase(sample_collection_date, report_date, infection.hospital());
            if (infection.fatal()) {
                long int medical_examiner_report_date = infection.deathTime + _par->deathReportingLag(REPORTING_RNG);
                // medical examiner only orders post-mortem PCR if sample was not already taken
                medical_examiner_report_date = max(medical_examiner_report_date, report_date);
                Community::reportDeath(sample_collection_date, medical_examiner_report_date);
            }
        }
    }
    // Flag locations with (non-historical) infections, so that we know to look there for human->mosquito transmission
    // Negative days are historical (pre-simulation) events, and thus we don't care about modeling transmission
    for (int day = std::max(infection.infectiousBegin, 0); day < infection.infectiousEnd; day++) {
        if (infection.inHospital(day)) {
            Location* hospital = getHospital();
            Community::flagInfectedLocation(this, infection.relInfectiousness, hospital->getType(), hospital, day);
        } else {
            Community::flagInfectedLocation(this, infection.relInfectiousness, getHomeLoc()->getType(), getHomeLoc(), day); // home loc can be a HOUSE or NURSINGHOME
            if (getDayLoc()) {
                // TODO -- people do not stop going to work/school when mild/moderately sick
                Community::flagInfectedLocation(this, infection.relInfectiousness, getDayLoc()->getType(), getDayLoc(), day);
            }
        }
    }


    // if the antibody-primed vaccine-induced immunity can be acquired retroactively, upgrade this person from naive to mature
    if (_par->retroactiveMatureVaccine) naiveVaccineProtection = false;

    return &infection;
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


bool Person::isCrossProtected(int time, StrainType strain) const { // assumes binary cross-immunity
    int cross_protected_until = INT_MIN;
    if (getNumNaturalInfections() > 0) {
        const Infection* inf = infectionHistory.back();
        //TODO -- Clean this up!
        const int base_days_immune = (signed) getDaysImmune();
        const int days_immune = inf->getStrain() == strain ? base_days_immune :
                                (inf->getStrain() == WILDTYPE and strain == B_1_1_7) or (inf->getStrain() == B_1_1_7 and strain == WILDTYPE) ? base_days_immune * 0.75 :
                                base_days_immune * 0.0; // one of the strains must be B_1_617_2
                                //(inf->getStrain() != B_1_617_2 and strain == B_1_617_2) or (inf->getStrain() == B_1_617_2 and strain != B_1_617_2) ? base_days_immune * 0.5;
        // cannot be reinfected until the last of {immunity, infectiousness, symptomatic period}
        vector<int> possible_end_dates = {inf->infectedBegin + days_immune, inf->infectiousEnd, inf->symptomEnd};
        cross_protected_until = covid::util::max_element(possible_end_dates);
    }
    return cross_protected_until >= time;
}


bool Person::isVaccineProtected(const int time, const StrainType strain) const {
    return isVaccinated() and
           ( !_par->vaccineLeaky or // if the vaccine isn't leaky
           (gsl_rng_uniform(RNG) < vaccineProtection(time, strain)) ); // or it protects (i.e., doesn't leak this time)
}


//bool Person::isSeroEligible(VaccineSeroConstraint vsc, double falsePos, double falseNeg) const {
bool Person::isSeroEligible() const {
    const VaccineSeroConstraint vsc = _par->vaccineSeroConstraint;
    const double falsePos = _par->seroTestFalsePos;
    const double falseNeg = _par->seroTestFalseNeg;

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
    if (isAlive(time)) {
        const size_t dose = vaccineHistory.size(); // this one isn't size() - 1, because it's the dose they're about to receive
        vaccineHistory.push_back(time);

        if ( isNaive() ) {
            naiveVaccineProtection = true;
        } else {
            naiveVaccineProtection = false;
        }

        if ( _par->vaccineLeaky == false ) { // all-or-none VE_S protection
            if ( (isNaive() and gsl_rng_uniform(RNG) < _par->VES_NAIVE_at(dose)) // vac someone who's naive
                 or gsl_rng_uniform(RNG) < _par->VES_at(dose) ) {                          // vac someone previously infected
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

void Person::selfQuarantine(const size_t today, const size_t quarantineDuration) {
    // set quarantine status
    quarantineStart = today;
    quarantineEnd   = today + quarantineDuration;
}

void Person::endQuarantine() {
    // clear quarantine status
    quarantineStart = INT_MAX;
    quarantineEnd   = INT_MAX;
}

bool Person::isQuarantining(const size_t today) {
    // check quarantine status
    return (quarantineStart != INT_MAX and quarantineEnd != INT_MAX) and (today >= quarantineStart) and (today <= quarantineEnd);
}
