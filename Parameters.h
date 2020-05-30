#ifndef __PARAMETERS_H
#define __PARAMETERS_H

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <bitset>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include "Utility.h"

enum TimePeriod {
    HOME,
    DAY,    // if this person leaves to go somewhere else, e.g. work or school
    NUM_OF_TIME_PERIODS
};

enum SexType {
    UNKNOWN,
    MALE,
    FEMALE,
    NUM_OF_SEX_TYPES
};

enum LocationType {
    HOUSE,
    WORK,
    SCHOOL,
    HOSPITAL,
    NURSINGHOME,
    NUM_OF_LOCATION_TYPES
};

enum ImmuneStateType {
    NAIVE,
    NATURAL,
    VACCINATED,
    NATURAL_AND_VACCINATED,
    NUM_OF_IMMUNE_STATE_TYPES
};

/*enum OutcomeType {
    ASYMPTOMATIC,
    MILD,
    SEVERE,
    CRITICAL,
    DEATH,
    NUM_OF_OUTCOME_TYPES
};*/

enum PathogenicityModel {
    CONSTANT_PATHOGENICITY,
    ORIGINAL_LOGISTIC,
    GEOMETRIC_PATHOGENICITY,
    NUM_OF_PRIMARY_PATHOGENICITY_MODELS
};

enum VaccineSeroConstraint {
    VACCINATE_SERONEGATIVE_ONLY,
    VACCINATE_SEROPOSITIVE_ONLY,
    VACCINATE_ALL_SERO_STATUSES,
    NUM_OF_VACCINE_SERO_CONSTRAINTS
};

extern const gsl_rng* RNG;// = gsl_rng_alloc (gsl_rng_taus2);

static const std::vector<std::string> MONTH_NAMES = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
static const std::vector<int> DAYS_IN_MONTH = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
static const std::vector<int> END_DAY_OF_MONTH = {31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};

static const int SYMPTOMATIC_DELAY = 1;                       // delay of symptoms after infectious period starts
//static const int INFECTIOUS_PERIOD_PRI = 5;                 // number of days until recovery from primary infection
//static const int INFECTIOUS_PERIOD_POST_PRI = 4;            // number of days until recovery from post-primary infection
static const int INFECTIOUS_PERIOD_ASYMPTOMATIC = 2;          // number of days until recovery for asymptomatic infections
static const int INFECTIOUS_PERIOD_MILD         = 4;          // number of days until recovery for mild cases
static const int INFECTIOUS_PERIOD_SEVERE       = 6;          // number of days until recovery for severe cases

// from Person.h
static const int NUM_AGE_CLASSES = 101;                       // maximum age+1 for a person
static const int MAX_INCUBATION = 9;                          // max incubation period in days
static const int MAX_HISTORY = 50;                            // length of exposure history in years

// cdf of incubation period, starting from day 1 (from Nishiura 2007)
const std::vector<double> INCUBATION_CDF = { 0,0,0.03590193,0.5070053,0.8248687,0.9124343,0.949212,0.974606,1 };

// for some serotypes, the fraction who are symptomatic upon primary infection
static const std::vector<double> SYMPTOMATIC_BY_AGE = {
    0.05189621,	0.05189621,	0.05189621,	0.05189621,	0.05189621,	0.1017964,	0.1017964,	0.1017964,   //  0-  7
    0.1017964,	0.1017964,	0.2774451,	0.2774451,	0.2774451,	0.2774451,	0.2774451,	0.4870259,   //  8- 15
    0.4870259,	0.4870259,	0.4870259,	0.4870259,	0.4870259,	0.4870259,	0.4870259,	0.4870259,   // 16- 23
    0.4870259,	0.8522954,	0.8522954,	0.8522954,	0.8522954,	0.8522954,	0.8522954,	0.8522954,   // 24- 31
    0.8522954,	0.8522954,	0.8522954,	0.9600798,	0.9600798,	0.9600798,	0.9600798,	0.9600798,   // 32- 39
    0.9600798,	0.9600798,	0.9600798,	0.9600798,	0.9600798,	1.0000000,	1.0000000,	1.0000000,   // 40- 47
    1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,   // 48- 55
    1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,   // 56- 59
    1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,   // 64- 71
    1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,   // 72- 79
    1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,   // 80- 87
    1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000,   // 88- 85
    1.0000000,	1.0000000,	1.0000000,	1.0000000,	1.0000000};                                      // 96-100

namespace covid {
    namespace standard {
        using std::cout;
        using std::cerr;
        using std::endl;
        using std::string;
        using std::vector;
        using std::set;
        using std::map;
        using std::pair;
        using std::make_pair;
        using std::ifstream;
        using std::istringstream;
        using std::ofstream;
        using std::ostringstream;
        using std::strcmp;
        using std::strtol;
        using std::strtod;
        using std::bitset;
    }
}


struct DynamicParameter {
    DynamicParameter(){};
    DynamicParameter(int s, int d, double v) : start(s), duration(d), value(v) {};
    int start;
    int duration;
    double value;
};


struct CatchupVaccinationEvent {
    CatchupVaccinationEvent(){};
    CatchupVaccinationEvent(int a, int s, double c): age(a), simDay(s), coverage(c) {};
    int age;
    int simDay;
    double coverage;
};


class Parameters {
public:

    Parameters() { define_defaults(); }
    Parameters(int argc, char *argv[]) { define_defaults(); readParameters(argc, argv); }

    void define_defaults();
    void readParameters(int argc, char *argv[]);
    void validate_parameters();
    void loadAnnualIntroductions(std::string annualIntrosFilename);

    static int sampler (const std::vector<double> CDF, const double rand, unsigned int index = 0) {
        while (index < CDF.size() and CDF[index] < rand) index++;
        return index;
    };

    unsigned long int randomseed;
    int runLength;
    double household_transmissibility;                      // per-day probability of transmission between a co-habitating dyad
    double workplace_transmissibility;                      // per-day probability of transmission, scaled by the fraction of infectious co-workers
    int infectiousOnset;                                    // days into infection when infectiousness starts
    int infectiousDuration;                                 // number of days infectious
    double VES;                                             // vaccine efficacy for susceptibility (can be leaky or all-or-none)
    double VES_NAIVE;                                       // VES for initially immunologically naive people
    double VEI;                                             // vaccine efficacy to reduce infectiousness
    double VEP;                                             // vaccine efficacy for pathogenicity
    double VEH;                                             // vaccine efficacy against hospitalization, given disease
    PathogenicityModel pathogenicityModel;                  // use age-specific values, or constant?
    double severeFraction;                                  // fraction of cases that are severe
    double criticalFraction;                                // fraction of severe cases that are critical
    double criticalMortality;                               // fraction of critical cases that die
    std::vector<double> hospitalizedFraction;               // Probability of being hospitalized, given asymptomatic, mild, severe, and critical infection
    std::vector<double> reportedFraction;                   // Probability of being reported, given asymptomatic, mild, severe, and critical infection
    bool vaccineLeaky;                                      // if false, vaccine is all-or-none
    bool retroactiveMatureVaccine;                          // if true, infection causes leaky vaccine to jump from naive to mature protection
    double seroTestFalsePos;                                // probability that seroneg person tests positive -- leaky test
    double seroTestFalseNeg;                                // probability that seropos person tests negative -- leaky test
    size_t numInitialExposed;                               // serotypes
    std::vector<double> numDailyExposed;                    // dimension is days
    size_t numInitialInfected;                              // serotypes
    double basePathogenicity;

    std::string populationFilename;
    std::string locationFilename;
    std::string peopleOutputFilename;
    std::string yearlyPeopleOutputFilename;
    std::string dailyOutputFilename;
    std::string annualIntroductionsFilename;                // time series of some external factor determining introduction rate
    std::vector<double> annualIntroductions;
    double annualIntroductionsCoef;                         // multiplier to rescale external introductions to something sensible
    int daysImmune;
    bool linearlyWaningVaccine;
    int vaccineImmunityDuration;
    bool vaccineBoosting;                                   // Are we re-vaccinated, either because of waning or because of multi-dose vaccine
    int numVaccineDoses;                                    // Number of times to boost; default is INT_MAX
    int vaccineDoseInterval;                                // How often to we re-vaccinate for initial vaccine course, in days
    int vaccineBoostingInterval;                            // How often to we re-vaccinate for boosting, in days
    std::vector<CatchupVaccinationEvent> catchupVaccinationEvents;
    int vaccineTargetAge;
    double vaccineTargetCoverage;
    int vaccineTargetStartDate;

    int startDayOfYear;
    bool dailyOutput;
    bool periodicOutput;
    int periodicOutputInterval;
    bool weeklyOutput;
    bool monthlyOutput;
    bool yearlyOutput;
    bool abcVerbose;
    unsigned long int serial;

    VaccineSeroConstraint vaccineSeroConstraint;
};

#endif
