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
    cerr << "\tis detected?     : " << boolalpha << (bool)_detection << noboolalpha        << endl;
    if (_detection) {
        cerr << "\t    detected state : " << _detection->detected_state                << endl;
        cerr << "\t    reported time  : " << _detection->reported_time                 << endl;
    }
}


Person::Person() {
    id = NEXT_ID++;
    age = -1;
    home_loc = nullptr;
    day_loc = nullptr;
    immune_state = NAIVE;
    //daysImmune = INT_MAX;
    naiveVaccineProtection = false;
    long_term_care = false;
    comorbidity = HEALTHY;
}


Person::~Person() {
    clearInfectionHistory();
}


void Person::revertState(const Date* date) {
    const int time = date->day();
    // remove infection history since date
    while (infectionHistory.size() > 0) {
        Infection* inf = infectionHistory.back();
        if (inf->getInfectedTime() >= time) {
            delete inf;
            infectionHistory.pop_back();
        } else {
            break;
        }
    }

    while (vaccineHistory.size() > 0) {
        if (vaccineHistory.back() >= time) {
            vaccineHistory.pop_back();
        } else {
            break;
        }
    }

    const bool had_infection   = infectionHistory.size() > 0;
    const bool had_vaccination = vaccineHistory.size() > 0;

    if (had_infection and had_vaccination) {
        immune_state = NATURAL_AND_VACCINATED;
    } else if (had_infection) {
        immune_state = NATURAL;
    } else if (had_vaccination) {
        immune_state = VACCINATED;
    } else {
        immune_state = NAIVE;
    }

    while (quarantineHistory.size() > 0) {
        if (quarantineHistory.back().first >= time) {
            quarantineHistory.pop_back();
        } else {
            break;
        }
    }
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


Infection& Person::initializeNewInfection(int time, size_t incubation_period, Location* sourceloc, Person* source) {
    Infection& infection      = initializeNewInfection();
    infection.infectionOwner  = this;
    infection.infectedBegin   = time;
    infection.infectedPlace   = sourceloc;
    infection.infectedBy      = source;
    infection.strain          = source ? source->getStrain() : WILDTYPE;
    infection.infectiousBegin = time + _par->infectiousness_onset(incubation_period);
    //infection.infectiousEnd   = time + incubation_period + SYMPTOMATIC_INFECTIOUS_PERIOD; // person may not actually be symptomatic!
    infection.infectiousEnd   = time + incubation_period + _par->strainPars[infection.strain].symptomaticInfectiousPeriod; // person may not actually be symptomatic!
//cerr << "iD " << infection.strain << ' ' << infection.infectiousBegin - time << ' ' << incubation_period << endl;
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


//bool Person::isInfectable(const int time, const StrainType strain) const {
//    return (not isInfected(time)) and                                   // not already infected
//           gsl_rng_uniform(RNG) < _par->susceptibilityByAge[age] and    // not innately resistant
//           !isCrossProtected(time, strain) and                          // no infection-based cross-immunity
//           !isVaccineProtected(time, strain);                           // no vaccine-based immunity
//}


//double Person::remainingEfficacy(const int time) const {
//    double remainingFraction = 1.0;
//    if (not isVaccinated()) {
//        remainingFraction = 0.0;
//    } else {
//        if (_par->linearlyWaningVaccine) {
//            // reduce by fraction of immunity duration that has waned
//            int time_since_vac = daysSinceVaccination(time);
//            if (time_since_vac > _par->vaccineImmunityDuration) {
//                remainingFraction = 0.0;
//            } else {
//                remainingFraction -= ((double) time_since_vac) / _par->vaccineImmunityDuration;
//            }
//        }
//    }
//    return remainingFraction;
//}


//// TODO -- merge this function with isVaccineProtected
//double Person::vaccineProtection(const int time, const StrainType strain) const {
//    double ves;
//    if (not isVaccinated()) {
//        ves = 0.0;
//    } else {
//        const size_t dose = vaccineHistory.size() - 1;
//        if (daysSinceVaccination(time) > _par->vaccineImmunityDuration) {
//            ves = 0.0;
//        } else {
//            if (naiveVaccineProtection == true) {
//                ves = _par->VES_NAIVE_at(dose, strain);
//            } else {
//                ves = _par->VES_at(dose, strain);
//            }
//            ves *= remainingEfficacy(time);
//        }
//    }
//    return ves;
//}

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


void Person::processDeath(Community* community, Infection &infection, const int deathTime) {
    community->tallyOutcome(DEATH);
    infection.deathTime         = deathTime;
    infection.infectiousEnd     = min(infection.infectiousEnd, deathTime);
    infection.symptomEnd        = min(infection.symptomEnd, deathTime);
    infection.severeEnd         = min(infection.severeEnd, deathTime);
    infection.criticalEnd       = min(infection.criticalEnd, deathTime);
}

// NOTE: PRNGs may be called (e.g., if immunity is leaky).  This function must be called for each unique exposure event.
vector<double> Person::calculate_event_probabilities(const int time, const StrainType strain) const {
    vector<double> Prob(NUM_OF_EVENT_PROBABILITY_TYPES, 0.0);
    const bool crossProtected   = isCrossProtected(time, strain);         // no infection-based cross-immunity
    const bool vaccineProtected = isVaccineProtected(time, strain);       // no vaccine-based immunity
    Prob[INFECTION_EVENT] = crossProtected or vaccineProtected ? 0.0 :_par->susceptibilityByAge[age];

    // current assumption is that only VES/IES wanes, other types of efficacy do not
    const size_t dose = vaccineHistory.size() - 1;
    const double effective_VEP = isVaccinated() ? _par->VEP_at(dose, strain) : 0.0;        // reduced pathogenicity due to vaccine
    const double effective_VEH = isVaccinated() ? _par->VEH_at(dose, strain) : 0.0;        // reduced severity (hospitalization) due to vaccine
    const double effective_VEF = isVaccinated() ? _par->VEF_at(dose, strain) : 0.0;        // reduced fatality due to vaccine

    double symptomatic_probability     = _par->pathogenicityByAge[age] * _par->strainPars[strain].relPathogenicity;
    symptomatic_probability           *= getNumNaturalInfections() ? (1.0 - _par->IEP) : 1.0;
    symptomatic_probability           *= 1.0 - effective_VEP;
    Prob[SYMPTOMATIC_EVENT] = symptomatic_probability;

    double severe_given_case           = _par->probSeriousOutcome.at(SEVERE)[comorbidity][age] * _par->strainPars[strain].relSeverity;
    severe_given_case                 *= getNumNaturalInfections() ? (1.0 - _par->IEH) : 1.0;
    severe_given_case                 *= 1.0 - effective_VEH;
    Prob[SEVERE_EVENT] = severe_given_case;

    const float hosp_prob = long_term_care ? LTC_SEVERE_TO_HOSPITAL : SEVERE_TO_HOSPITAL;
    Prob[HOSPITALIZATION_EVENT] = hosp_prob;

    const double critical_given_severe = _par->probSeriousOutcome.at(CRITICAL)[comorbidity][age];
    Prob[CRITICAL_EVENT] = critical_given_severe;

    double mortalityCoef               = _par->strainPars[strain].relMortality;
    mortalityCoef                     *= getNumNaturalInfections() ? (1.0 - _par->IEF) : 1.0;
    mortalityCoef                     *= 1.0 - effective_VEF;

    const double nonIcuMotality        = NON_ICU_CRITICAL_MORTALITY * mortalityCoef;        // icu mortality calculated later, as it depends on timing
    Prob[NON_ICU_DEATH_EVENT] = nonIcuMotality;

    Prob[ICU_DEATH_EVENT_COEF] = _par->strainPars[strain].relIcuMortality * mortalityCoef;
    return Prob;
}

// infect - infect this individual
// returns non-null pointer if infection occurs
Infection* Person::infect(Community* community, Person* source, const Date* date, Location* sourceloc, StrainType strain, bool /*check_susceptibility*/) {
    const int time = date->day();
    const vector<double> Pr = calculate_event_probabilities(time, strain);
    // Bail now if this person can not become infected
    // Not quite the same as "susceptible"--this person may be e.g. partially immune
    // due to natural infection or vaccination
    if (isInfected(time) or gsl_rng_uniform(RNG) > Pr[INFECTION_EVENT]) { return nullptr; }

    // Create a new infection record
    const size_t incubation_period = _par->symptom_onset(strain); // may not be symptomatic, but this is used to determine infectiousness onset
    Infection& infection = initializeNewInfection(time, incubation_period, sourceloc, source);
    community->tallyOutcome(ASYMPTOMATIC);
    if (not source) { infection.strain = strain; }

    const size_t dose = vaccineHistory.size() - 1;
    const double effective_VEI = isVaccinated() ? _par->VEI_at(dose, strain) : 0.0;        // reduced infectiousness
    infection.relInfectiousness       *= _par->strainPars[strain].relInfectiousness;
    infection.relInfectiousness       *= getNumNaturalInfections() > 1 ? (1.0 - _par->IEI) : 1.0; // getNumNaturalInfections() counts this infection too
    infection.relInfectiousness       *= 1.0 - effective_VEI;

    const double highly_infectious_threshold = 8.04; // 80th %ile for overall SARS-CoV-2 from doi: 10.7554/eLife.65774, "Fig 4-Fig Sup 3"
    const double rel_infectiousness_high     = 4.0;
    const double rel_infectiousness_norm     = 0.25;

    const double asymp_weibull_scale = 6.72;
    const double asymp_weibull_exp   = 3.33;
    const double symp_weibull_scale  = 7.40;
    const double symp_weibull_exp    = 3.81;

    // determine disease outcome and timings
    if ( gsl_rng_uniform(RNG) > Pr[SYMPTOMATIC_EVENT]) {
       // asymptomatic; est mean relInfectiousness = 6.03
        infection.relInfectiousness *= gsl_ran_weibull(RNG, asymp_weibull_scale, asymp_weibull_exp) > highly_infectious_threshold ? rel_infectiousness_high : rel_infectiousness_norm;
    } else {
        // symptomatic; est mean relInfectiousness = 6.69
        infection.relInfectiousness *= gsl_ran_weibull(RNG, symp_weibull_scale, symp_weibull_exp) > highly_infectious_threshold ? rel_infectiousness_high : rel_infectiousness_norm;
        //const size_t symptom_onset = _par->symptom_onset();
        community->tallyOutcome(MILD);
        infection.symptomBegin = time + incubation_period;
        if ( not (gsl_rng_uniform(RNG) < Pr[SEVERE_EVENT]) ) {
            // It does not become severe
            infection.symptomEnd = infection.symptomBegin + _par->symptom_duration_mild();
        } else {
            // It does progress and become severe
            community->tallyOutcome(SEVERE);
            infection.severeBegin = infection.symptomBegin + _par->pre_severe_symptomatic();

            // Is this person hospitalized when their severe symptoms begin?
            bool hosp = false;
            if (gsl_rng_uniform(RNG) < Pr[HOSPITALIZATION_EVENT]) {
                // This person is hospitalized when symptoms become severe
                hosp = true;
                infection.hospitalizedBegin = infection.severeBegin;
            }

            if (not (gsl_rng_uniform(RNG) < Pr[CRITICAL_EVENT])) {
                // Severe, but does not become critical
                infection.severeEnd     = infection.severeBegin   + _par->severe_only_duration();
                infection.symptomEnd    = infection.severeEnd; // TODO - extend symptoms beyond severe period (relevant for e.g. econ analyses)
            } else {
                // It does progress to critical disease
                community->tallyOutcome(CRITICAL);
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
                    const double icuMortality          = _par->icuMortality(comorbidity, age, infection.icuBegin) * Pr[ICU_DEATH_EVENT_COEF];

                    if (not hosp) { infection.hospitalizedBegin = infection.icuBegin; } // if they weren't hospitalized before, they are now
                    death = gsl_rng_uniform(RNG) < icuMortality;
                    if (death) {
                        // uniform randomly chose a day from the critical duration when death happens
                        processDeath(community, infection, infection.criticalBegin + _par->sampleIcuTimeToDeath());
                    }
                } else {
                    death = gsl_rng_uniform(RNG) < Pr[NON_ICU_DEATH_EVENT];
                    if (death) {
                        // non-icu death, happens when critical symptoms begin, as this person is not receiving care
                        processDeath(community, infection, infection.criticalBegin + _par->sampleCommunityTimeToDeath());
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

        if (report_date >= (long int) _par->runLength) { detected = false; }
        // this needs to be evaluated early to ensure consistent number of REPORTING_RNG draws regardless of runLength
        size_t death_rep_lag = infection.fatal() ? _par->deathReportingLag(REPORTING_RNG) : 0;

        if (detected) {
            infection.detect(detected_state, report_date);
            community->reportCase(sample_collection_date, report_date, infection.hospital());
            if (infection.fatal()) {
                long int medical_examiner_report_date = infection.deathTime + death_rep_lag;
                // medical examiner only orders post-mortem PCR if sample was not already taken
                medical_examiner_report_date = max(medical_examiner_report_date, report_date);
                //community->reportDeath(sample_collection_date, medical_examiner_report_date);
                community->reportDeath(infection.deathTime, medical_examiner_report_date);
            }
        }
    }
    // Flag locations with (non-historical) infections, so that we know to look there for human->mosquito transmission
    // Negative days are historical (pre-simulation) events, and thus we don't care about modeling transmission
    for (int day = std::max(infection.infectiousBegin, 0); day < infection.infectiousEnd; day++) {
        if (infection.inHospital(day)) {
            Location* hospital = getHospital();
            community->flagInfectedLocation(this, infection.relInfectiousness, hospital->getType(), hospital, day);
        } else {
            community->flagInfectedLocation(this, infection.relInfectiousness, getHomeLoc()->getType(), getHomeLoc(), day); // home loc can be a HOUSE or NURSINGHOME
            if (getDayLoc() and not (infection.isSevere(day) or infection.isCritical(day))) {
                // TODO -- people do not stop going to work/school when mild/moderately sick
                community->flagInfectedLocation(this, infection.relInfectiousness, getDayLoc()->getType(), getDayLoc(), day);
            }
        }
    }

    // community->tallyInfectionByLoc(&infection);

    // if the antibody-primed vaccine-induced immunity can be acquired retroactively, upgrade this person from naive to mature
    if (_par->retroactiveMatureVaccine) naiveVaccineProtection = false;

    return &infection;
}


bool Person::isCrossProtected(int time, StrainType strain) const { // assumes binary cross-immunity
    bool immune = false;

    if (getNumNaturalInfections() > 0) {
        // first check for broad, short term immunity
        const Infection* last_inf = infectionHistory.back();
        vector<int> possible_last_infection_end_dates = {last_inf->infectiousEnd, last_inf->symptomEnd};
        const int last_infection_end_date = covid::util::max_element(possible_last_infection_end_dates);
        const int time_since_last_infection = time - last_infection_end_date;
        double remaining_natural_efficacy = _par->remainingEfficacy(startingNaturalEfficacy, time_since_last_infection);
        remaining_natural_efficacy *= 1.0 - _par->strainPars[strain].immuneEscapeProb;
        if ( (_par->immunityLeaky and (gsl_rng_uniform(RNG) < remaining_natural_efficacy)) or
             (not _par->immunityLeaky and (time_since_last_infection < naturalImmunityDuration)) ) {
            immune = true;
        }

        if (not immune) {
            // maybe there's longer term, strain-specific immunity
            for (Infection* inf: infectionHistory) {
                const int strain_sp_immune_duration     = 60;
                const double strain_sp_immune_leakiness = 0.85;
                //if (inf->getStrain() == strain and
                if (_par->crossProtectionMatrix[inf->getStrain()][strain] and
                    (time_since_last_infection < strain_sp_immune_duration or gsl_rng_uniform(RNG) < strain_sp_immune_leakiness)) { // CABP: https://www.sciencedirect.com/science/article/pii/S0140673621006759
                    immune = true;
                    break;
                }
            }
        }

        if (not immune) {
            // longer-term broad immunity
            if (crossProtectionProbability < _par->strainPars[strain].immuneEscapeProb) {
                immune = true;
            }
        }
    }
    return immune;
}


bool Person::isVaccineProtected(const int time, const StrainType strain) const {
    bool immune = false;
    if (isVaccinated()) {
        const size_t dose = vaccineHistory.size() - 1;
        double init_ves;
        if (naiveVaccineProtection == true) {
            init_ves = _par->VES_NAIVE_at(dose, strain);
        } else {
            init_ves = _par->VES_at(dose, strain);
        }

        double remaining_vaccine_efficacy = _par->remainingEfficacy(init_ves, daysSinceVaccination(time) - _par->vaccine_dose_to_protection_lag);
        remaining_vaccine_efficacy *= 1.0 - _par->strainPars[strain].immuneEscapeProb;

        const bool beyond_vaccine_dose_protection_lag = daysSinceVaccination(time) >= _par->vaccine_dose_to_protection_lag;
        const bool protection_remains = _par->immunityLeaky ? gsl_rng_uniform(RNG) < remaining_vaccine_efficacy
                                                           : (daysSinceVaccination(time) - _par->vaccine_dose_to_protection_lag) < getVaccineImmunityDuration(dose, strain);

        immune = beyond_vaccine_dose_protection_lag and protection_remains;
    }
    return immune;
}


bool Person::isSeroEligible() const {
    const VaccineSeroConstraint vsc = _par->vaccineSeroConstraint;
    const double falsePos = _par->seroTestFalsePos;
    const double falseNeg = _par->seroTestFalseNeg;

    if (vsc == VACCINATE_ALL_SERO_STATUSES) return true;

    assert(falsePos >= 0.0 and falsePos <= 1.0);
    assert(falseNeg >= 0.0 and falseNeg <= 1.0);
    assert(vsc == VACCINATE_SEROPOSITIVE_ONLY or vsc == VACCINATE_SERONEGATIVE_ONLY);

    bool isSeroPos = not isImmuneState(NAIVE); // fully susceptible == seronegative == false

    if ((isSeroPos and (falseNeg > gsl_rng_uniform(RNG)))       // sero+ but tests negative
        or (!isSeroPos and (falsePos > gsl_rng_uniform(RNG)))) { // sero- but tests positive
        isSeroPos = !isSeroPos;
    }

    bool eligible = vsc == VACCINATE_SEROPOSITIVE_ONLY ? isSeroPos : not isSeroPos;

    return eligible;
}


bool Person::isInfEligible(int today) const {
    const VaccineInfConstraint vic = _par->vaccineInfConstraint;

    if (vic == VACCINATE_ALL_INF_STATUSES) return true;

    const bool prior_inf = infectionHistory.size() ? true : false;
    bool prior_case = false;
    for (Infection* inf : infectionHistory) {
        if (inf->isDetected() and inf->getDetection()->reported_time <= today) { prior_case = true; break; }
    }

    bool eligible = true;

    switch (vic) {
        case VACCINATE_NAIVE_ONLY:       { eligible = not prior_inf; break; } // no prior infections
        case VACCINATE_INF_ONLY:         { eligible = prior_inf; break; }  // any prior infections
        case VACCINATE_NON_CASE_ONLY:    { eligible = not prior_case; break;} // no prior detected infections
        case VACCINATE_CASE_ONLY:        { eligible = prior_case; break; } // any prior detected infections
        default:                         { /*should never be here*/ break; }
    }

    return eligible;
}


bool Person::vaccinate(int time) {
    if (isAlive(time)) {
        const size_t dose = vaccineHistory.size(); // this one isn't size() - 1, because it's the dose they're about to receive
        vaccineHistory.push_back(time);

        if ( isImmuneState(NAIVE) ) {
            naiveVaccineProtection = true;
        } else {
            naiveVaccineProtection = false;
        }

        assert( _par->immunityLeaky == true ); // Currently we do not support all-or-none immunity, b/c of how we model reduced protection due to VOCs
        if ( _par->immunityLeaky == true or      // vaccine is leaky, so update immune state
           ( _par->immunityLeaky == false and    // vaccine is all-or-none, so update immune state if...
               ( (isImmuneState(NAIVE) and gsl_rng_uniform(RNG) < _par->VES_NAIVE_at(dose))             // person is naive and naive VES takes
                   or (not isImmuneState(NAIVE) and gsl_rng_uniform(RNG) < _par->VES_at(dose)) ) ) ) {  // person is not naive and non-naive VES takes

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
        return true;
    } else {
        return false;
    }
}

void Person::scheduleQuarantine(int time, int quarantineDuration) {
    quarantineHistory.emplace_back(time, time + quarantineDuration);
}

bool Person::isQuarantining(int time) {
    if (quarantineHistory.size() > 0) {
        const pair<int, int> interval = quarantineHistory.back();
        return time >= interval.first and time < interval.second;
    } else {
        return false;
    }
}
