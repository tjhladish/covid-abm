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
#include <climits>
#include "Utility.h"

class Date;

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

enum OutcomeType {
    ASYMPTOMATIC,
    MILD,
    SEVERE,
    CRITICAL,
    DEATH,
    NUM_OF_OUTCOME_TYPES
};

enum ComorbidType{
    HEALTHY,
    COMORBID,
    NUM_OF_COMORBID_TYPES
};

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

enum TimedIntervention {
    SCHOOL_CLOSURE,
    NONESSENTIAL_BUSINESS_CLOSURE,
    SOCIAL_DISTANCING,
    NUM_OF_TIMED_INTERVNETIONS
};

enum MmodsScenario {
    MMODS_CLOSED,   //0
    MMODS_2WEEKS,   //1
    MMODS_5PERCENT, //2
    MMODS_OPEN,     //3
    NUM_OF_MMODS_SCENARIOS
};

struct GammaPars {
    GammaPars() {};
    GammaPars(double _a, double _b) : a(_a), b(_b) {};
    // gamma distribution parameterization used by GSL
    double a; // shape
    double b; // scale
};

class ReportingLagModel {
  public:
    ReportingLagModel() {};
    ReportingLagModel(std::string filename) { read_csv(filename); };
    void insert_lag (std::string ymd, double a_shape, double b_scale) {
        dated_lags[ymd] = GammaPars(a_shape, b_scale);
    }
    std::string earliest_date() const { return dated_lags.begin()->first; }
    std::string latest_date()   const { return dated_lags.rbegin()->first; }

    void read_csv(std::string filename) {
        std::vector<std::vector<std::string>> reporting_lag_data = covid::util::read_2D_vector_file(filename, ',');
        bool header = true;
        for (size_t i = header; i < reporting_lag_data.size(); ++i) {
            vector<string> row = reporting_lag_data[i];
            assert(row.size() == 5);
            string date    = row[0];
            double a_shape = stod(row[3]);
            double b_scale = stod(row[4]);
            insert_lag(date, a_shape, b_scale);
        }

    }

    size_t sample(const gsl_rng* REPORTING_RNG, std::string ymd) const {
        if (ymd < earliest_date()) {
            return INT_MAX;
        } else if (ymd > latest_date()) {
            GammaPars gp = dated_lags.rbegin()->second;
            return (size_t) round(gsl_ran_gamma(REPORTING_RNG, gp.a, gp.b));
        } else if (dated_lags.count(ymd) > 0) {
            GammaPars gp = dated_lags.at(ymd);
            return (size_t) round(gsl_ran_gamma(REPORTING_RNG, gp.a, gp.b));
        } else {
            cerr << "ERROR: The date provided to ReportingLagModel::sample() (" << ymd
                 << ") is in the range of known dates, but lag parameters are not known for that date.\n";
            exit(-1);
        }
    }

    size_t sample(const gsl_rng* REPORTING_RNG, const Date* date) const;

  private:
    std::map<std::string, GammaPars> dated_lags;
};


extern const gsl_rng* RNG;// = gsl_rng_alloc (gsl_rng_taus2);
// use a second RNG for stochastic reporting; this allows for more powerful analysis of
// the effects of different reporting models
extern const gsl_rng* REPORTING_RNG; // = gsl_rng_alloc (gsl_rng_mt19937);

/*
// transmission-related probabilities
susceptibility -- use fuction wrapper for generality so we can potentially use age in the future, but just a single value in par for now
can be 1.0 for right now

transmission -- wrapper function for setting and age, with covariates of 1 for now
trans per hour probability to be fit
hh secondary infection probability between 15-30%
transmission from asymp is 0.5*trans from symp

// disease-related probabilities
/// probability of outcomes
-- infection -- to be fit/prior informed by hh secondary attack rate
-- symptomatic -- roughly 0.5, but should be age-structured (and we have some estimates)
-- severe (hosp is a subset) -- roughly 0.1 (severe given symptoms), also should be age-structured
    - fraction of severe that are hospitalized = 0.8 general pop, 0.5 of nursinghome pop
-- critical (icu is a subset) -- 0.5 (critical given severe)
    - critical-->icu if hospitalized, 1.0; 0.5 otherwise
-- death -- IFR ~ 0.005, 0.2 of critical infections in icu result in death, 1.0 of non-icu critical

/// times
.time to infectiousness ~3 days
.duration of infectiousness 7 days
.time to symptoms ~5 days from infection
.time to severe ~7 days from symptoms
.time to critical 4 days from severe
--the previous numbers regarding severe disease surrounding critical disease contradict the following
.duration of symptoms 14 days for mild, (1 week mild, 2w severe, 1w mild) for severe, 10 more days critical, in the middle of severe period
*/

//static const float SEVERE_FRACTION = 0.1;                   // fraction of cases become severe
//static const float CRITICAL_FRACTION = 0.33;                  // fraction of severe cases that become critical
                                                              // https://www.cdc.gov/mmwr/volumes/69/wr/mm6932e3.htm?s

static const float SEVERE_TO_HOSPITAL = 0.8;                  // general population, probability of going to hospital if severe
static const float LTC_SEVERE_TO_HOSPITAL = 0.25;             // probability long-term-care residents who are severe go to hospital (if hosp available)
static const float CRITICAL_TO_ICU_IF_HOSP = 0.9;             // probability already-hospitalized patients go to ICU when critical (if ICU available)
static const float CRITICAL_TO_ICU_IF_NOT_HOSP = 0.75;         // probability non-hospitalized severe patients go to ICU when critical (if ICU available)

static const int INFECTIOUSNESS_ONSET = 3;                    // delay of infectiousness after infectious exposure
static const int INFECTIOUS_PERIOD = 7;                       // independent of disease outcome
static const int SYMPTOM_ONSET = 5;                           // delay of symptoms after infectious exposure
static const int SYMPTOM_DURATION_MILD = 14;                  // number of days symptoms remain for people who do not become severe

// severe symptoms warrant hospitalization (which may or may not occur)
static const int PRE_SEVERE_SYMPTOMATIC = 7;                  // when severe disease occurs, delay from symptom onset
static const int SEVERE_DURATION = 14;                        // duration of severe disease
static const int POST_SEVERE_SYMPTOMATIC = 7;                 // number of days after severe disease that symptoms remain
static const int SYMPTOM_DURATION_SEVERE = PRE_SEVERE_SYMPTOMATIC + SEVERE_DURATION + POST_SEVERE_SYMPTOMATIC;

// critical symptoms warrant intensive care (which may or may not occur)
static const int PRE_CRITICAL_SEVERE = 4;                     // when critical disease occurs, delay from severe onset
static const int CRITICAL_DURATION = 10;                      // duration of critical disease

static const int POST_CRITICAL_SEVERE = 10;                   // number of days after critical disease that severe symptoms remain
static const int SEVERE_DURATION_CRITICAL = PRE_CRITICAL_SEVERE + CRITICAL_DURATION + POST_CRITICAL_SEVERE;
static const int SYMPTOM_DURATION_CRITICAL = PRE_SEVERE_SYMPTOMATIC + SEVERE_DURATION_CRITICAL + POST_SEVERE_SYMPTOMATIC;
//static const float ICU_CRITICAL_MORTALITY = 0.4;              // baseline probability of dying while in ICU, distributed uniformly across days in ICU
                                                              // 0.4 based on ARDS from here https://ccforum.biomedcentral.com/articles/10.1186/s13054-020-03006-1
static const float NON_ICU_CRITICAL_MORTALITY = 1.0;          // probability of dying when becoming critical if intensive care is not received

// from Person.h
static const int NUM_AGE_CLASSES = 121;                       // maximum age+1 for a person


// for some serotypes, the fraction who are symptomatic upon primary infection
/*static const std::vector<double> SYMPTOMATIC_BY_AGE = {
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
*/


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
    ~Parameters() { if (rlm) delete rlm; }

    void define_defaults();
    void readParameters(int argc, char *argv[]);
    void validate_parameters();
    void loadAnnualIntroductions(std::string annualIntrosFilename);
    void define_susceptibility_and_pathogenicity();

    static int sampler (const std::vector<double> CDF, const double rand, unsigned int index = 0) {
        while (index < CDF.size() and CDF[index] < rand) index++;
        return index;
    };

    static vector<double> toReportedFraction(vector<double> rho_detect) {
        vector<double> reported_fraction;
        for (double rho: rho_detect) {
            if (reported_fraction.size() == 0) {
                reported_fraction.push_back(rho);
            } else {
                reported_fraction.push_back(1.0 - (1.0 - rho)*(1.0 - reported_fraction.back()));
            }
        }
        return reported_fraction;
    }

    int sampleIcuTimeToDeath() const {
    // parameterization based on what you need to produce findings of median = 7, IQR = [3,11]
    // https://www.thelancet.com/journals/lanres/article/PIIS2213-2600(20)30079-5/fulltext
        const double prob_of_daily_survival = 0.2;
        const size_t num_days_survive_until_death = 2;
        return gsl_ran_negative_binomial(RNG, prob_of_daily_survival, num_days_survive_until_death);
    }

    double timedInterventionEffect(TimedIntervention ti, size_t day) const;

    unsigned long int randomseed;
    size_t runLength;
    double household_transmissibility;                      // per-day probability of transmission between a co-habitating dyad
    double workplace_transmissibility;                      // per-day probability of transmission, scaled by the fraction of infectious co-workers
    double school_transmissibility;                         // per-day probability of transmission, scaled by the fraction of infectious students/staff
    double social_transmissibility;                         // per-day probability of transmission, scaled by the fraction of infectious social contacts
    vector<float> susceptibilityByAge;                      // probability of infection given exposure, index by year of age
    vector<float> pathogenicityByAge;                       // probability of clinical disease given infection, index by year of age
    vector<float> severeFractionByAge;                      // probability of severe disease given clinical disease, index by year of age
    map<OutcomeType, vector< vector<float>>> probSeriousOutcome;        // look up probability of severe, critical, fatal outcomes by age & comorbid status
    double VES;                                             // vaccine efficacy for susceptibility (can be leaky or all-or-none)
    double VES_NAIVE;                                       // VES for initially immunologically naive people
    double VEI;                                             // vaccine efficacy to reduce infectiousness
    double VEP;                                             // vaccine efficacy for pathogenicity
    double VEH;                                             // vaccine efficacy against hospitalization, given disease
    PathogenicityModel pathogenicityModel;                  // use age-specific values, or constant?
    //std::vector<double> reportedFraction;                   // Prob of being reported, given asymptomatic, mild, severe, critical, and fatal infection
    std::vector<double> probFirstDetection;                 // Probability of being *first* detected while {asymp, mild, severe, crit, dead}
    size_t symptomToTestLag;                                // For people with mild disease who get infected, number of days until tested
    void createReportingLagModel(std::string filename);
    size_t defaultReportingLag;                             // Number of days it takes for a test result to be reported after a sample is collected;
                                                            // only used if an explicit reporting lag model is not defined
    ReportingLagModel* rlm;
    size_t reportingLag (const gsl_rng* REPORTING_RNG, const Date *date) const {
        if (rlm) {
            return rlm->sample(REPORTING_RNG, date);
        } else {
            return defaultReportingLag;
        }
    }

    void createIcuMortalityReductionModel(double maximum_val, double inflection_sim_day, double slope);
    double icuMortality(ComorbidType comorbidity, size_t age, size_t sim_day) const;
    double icuMortalityFraction;                            // fraction of all deaths that occur in ICUs;
                                                            // used for interpreting empirical mortality data, *not within simulation*

    size_t deathReportingLag;                               // number of days from when death occurs to when it's reported
    void createSocialDistancingModel(std::string filename, float mobility_logit_shift, float mobility_logit_stretch);
    bool vaccineLeaky;                                      // if false, vaccine is all-or-none
    bool retroactiveMatureVaccine;                          // if true, infection causes leaky vaccine to jump from naive to mature protection
    double seroTestFalsePos;                                // probability that seroneg person tests positive -- leaky test
    double seroTestFalseNeg;                                // probability that seropos person tests negative -- leaky test
    size_t numInitialExposed;                               // serotypes
//    std::vector<double> numDailyExposed;                    // dimension is days
    std::vector<double> probDailyExposure;                  // per person, per day probability of exposure
    size_t numInitialInfected;                              // serotypes

    std::string populationFilename;
    std::string comorbidityFilename;
    std::string locationFilename;
    std::string networkFilename;
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
    size_t numSurveilledPeople;

    bool traceContacts;                                     // should we keep track of who infected whom

                                                            // e.g. school closures, non-essential business closures, social distancing
    std::map<TimedIntervention, std::vector<float>> timedInterventions;
    std::vector<double> icuMortalityReduction;              // time-varying reduction (from 0) in ICU mortality due to improved Tx

    size_t startDayOfYear;
    size_t julianYear;
    bool dailyOutput;
    bool periodicOutput;
    int periodicOutputInterval;
    bool weeklyOutput;
    bool monthlyOutput;
    bool yearlyOutput;
    bool abcVerbose;
    unsigned long int serial;

    VaccineSeroConstraint vaccineSeroConstraint;
    MmodsScenario mmodsScenario;
};

#endif
