#ifndef __PARAMETERS_H
#define __PARAMETERS_H
#define varname(s) #s

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
#include <gsl/gsl_cdf.h>
#include <climits>
#include "Utility.h"

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

class Date;

extern gsl_rng* RNG;// = gsl_rng_alloc (gsl_rng_taus2);
// use a second RNG for stochastic reporting; this allows for more powerful analysis of
// the effects of different reporting models
extern gsl_rng* REPORTING_RNG; // = gsl_rng_alloc (gsl_rng_mt19937);

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

inline std::ostream& operator<<(std::ostream& out, const SexType value){
    const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value){
        PROCESS_VAL(UNKNOWN);
        PROCESS_VAL(MALE);
        PROCESS_VAL(FEMALE);
        PROCESS_VAL(NUM_OF_SEX_TYPES);
    }
#undef PROCESS_VAL
    return out << s;
}

enum LocationType {
    HOUSE,
    WORK,
    SCHOOL,
    HOSPITAL,
    NURSINGHOME,
    NUM_OF_LOCATION_TYPES
};

inline std::ostream& operator<<(std::ostream& out, const LocationType value){
    const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value){
        PROCESS_VAL(HOUSE);
        PROCESS_VAL(WORK);
        PROCESS_VAL(SCHOOL);
        PROCESS_VAL(HOSPITAL);
        PROCESS_VAL(NURSINGHOME);
        PROCESS_VAL(NUM_OF_LOCATION_TYPES);
    }
#undef PROCESS_VAL
    return out << s;
}

enum PublicTransmissionType {
    NO_PUBLIC_TRANSMISSION,
    LOW_PUBLIC_TRANSMISSION,
    HIGH_PUBLIC_TRANSMISSION,
    NUM_OF_PUBLIC_TRANSMISSION_TYPES
};

inline std::ostream& operator<<(std::ostream& out, const PublicTransmissionType value){
    const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value){
        PROCESS_VAL(NO_PUBLIC_TRANSMISSION);
        PROCESS_VAL(LOW_PUBLIC_TRANSMISSION);
        PROCESS_VAL(HIGH_PUBLIC_TRANSMISSION);
        PROCESS_VAL(NUM_OF_PUBLIC_TRANSMISSION_TYPES);
    }
#undef PROCESS_VAL
    return out << s;
}

// immune state does not indicate protection, only history of infection/vaccination
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

inline std::ostream& operator<<(std::ostream& out, const OutcomeType value){
    const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value){
        PROCESS_VAL(ASYMPTOMATIC);
        PROCESS_VAL(MILD);
        PROCESS_VAL(SEVERE);
        PROCESS_VAL(CRITICAL);
        PROCESS_VAL(DEATH);
        PROCESS_VAL(NUM_OF_OUTCOME_TYPES);
    }
#undef PROCESS_VAL
    return out << s;
}

enum ComorbidType{
    HEALTHY,
    COMORBID,
    NUM_OF_COMORBID_TYPES
};

enum StrainType{
    WILDTYPE,
    ALPHA, //B_1_1_7,
    DELTA, //B_1_617_2,
    OMICRON,
    NUM_OF_STRAIN_TYPES
};

inline std::ostream& operator<<(std::ostream& out, const StrainType value){
    const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value){
        PROCESS_VAL(WILDTYPE);
        PROCESS_VAL(ALPHA);
        PROCESS_VAL(DELTA);
        PROCESS_VAL(OMICRON);
        PROCESS_VAL(NUM_OF_STRAIN_TYPES);
    }
#undef PROCESS_VAL
    return out << s;
}

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

enum CsmhScenario {
    CSMH_A,
    CSMH_B,
    CSMH_C,
    CSMH_D,
    NUM_OF_CSMH_SCENARIOS
};

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

static const float SEVERE_TO_HOSPITAL = 0.52; // tried 0.1     // general population, probability of going to hospital if severe
static const float LTC_SEVERE_TO_HOSPITAL = 0.2; // 0.1       // probability long-term-care residents who are severe go to hospital (if hosp available)
static const float CRITICAL_TO_ICU_IF_HOSP = 1.0;             // probability already-hospitalized patients go to ICU when critical (if ICU available)
static const float CRITICAL_TO_ICU_IF_NOT_HOSP = 0.5;         // probability non-hospitalized severe patients go to ICU when critical (if ICU available)

// symptom onset gamma parameters: https://elifesciences.org/articles/57149
static const float WILDTYPE_SYMPTOM_ONSET_MEAN = 4.91;
static const float WILDTYPE_SYMPTOM_ONSET_GAMMA_SHAPE = 3.05;                                                               // shape parameter for sampling incubation period
static const float WILDTYPE_SYMPTOM_ONSET_GAMMA_SCALE = WILDTYPE_SYMPTOM_ONSET_MEAN/WILDTYPE_SYMPTOM_ONSET_GAMMA_SHAPE;     // scale parameter for sampling incubation period

// number of days symptoms remain for people who do not become severe
// fit to data from https://doi.org/10.1111/joim.13089
static const float WILDTYPE_SYMPTOM_DURATION_MILD_GAMMA_SHAPE = 4.070483;
static const float WILDTYPE_SYMPTOM_DURATION_MILD_GAMMA_SCALE = 2.825217;

// expected fraction of incubation period spent non-infectious, used to sample from binomial
// https://www.ams.edu.sg/view-pdf.aspx?file=media%5c5556_fi_331.pdf&ofile=Period+of+Infectivity+Position+Statement+(final)+23-5-20+(logos).pdf
// (↑↑↑ ref for 2.3 day asymptomatic transmission mean)
static const float ASYMPTOMATIC_TRANSMISSION_MEAN = 2.3;
static const float INFECTIOUSNESS_ONSET_FRACTION  = 1.0 - (ASYMPTOMATIC_TRANSMISSION_MEAN/WILDTYPE_SYMPTOM_ONSET_MEAN);

// https://www.ams.edu.sg/view-pdf.aspx?file=media%5c5556_fi_331.pdf&ofile=Period+of+Infectivity+Position+Statement+(final)+23-5-20+(logos).pdf
static const int SYMPTOMATIC_INFECTIOUS_PERIOD = 7;           // independent of disease outcome

// severe symptoms warrant hospitalization (which may or may not occur)
// when severe disease occurs, delay from symptom onset
// fit to data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7171478/
static const float PRE_SEVERE_SYMPTOMATIC_GAMMA_SHAPE = 0.812130;
static const float PRE_SEVERE_SYMPTOMATIC_GAMMA_SCALE = 7.199686;

// Empirical hospital lenght-of-stay calculated based on github data (distribution_general_world.csv, distribution_icu_world.csv) provided by
// https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01726-3#Sec13
// As that data was for total hospital stay, we subtracted off the density for ICU stays, assuming that 1/3 of hospital stays result in ICU
// stays (https://www.cdc.gov/mmwr/volumes/69/wr/mm6932e3.htm?s).
static const vector<double> SEVERE_ONLY_DURATION_CDF =
    {0.037435, 0.139665, 0.256730, 0.374955, 0.491195, 0.594245, 0.683795, 0.757310, 0.812225, 0.855110,
     0.886960, 0.909490, 0.925635, 0.935870, 0.943265, 0.949685, 0.954360, 0.958955, 0.963360, 0.967855,
     0.971540, 0.975125, 0.978495, 0.981425, 0.984275, 0.986235, 0.988115, 0.989515, 0.991115, 0.992335,
     0.993585, 0.994790, 0.995595, 0.996245, 0.996965, 0.997490, 0.997980, 0.998385, 0.998670, 0.998810,
     0.999055, 0.999170, 0.999255, 0.999305, 0.999335, 0.999420, 0.999545, 0.999585, 0.999585, 0.999645,
     0.999660, 0.999740, 0.999830, 0.999905, 0.999920, 0.999935, 0.999960, 0.999970, 0.999980, 0.999985,
     0.999985, 0.999985, 0.999985, 0.999985, 0.999985, 0.999985, 0.999985, 0.999985, 0.999985, 0.999985,
     1.000000};

// As above, empirical ICU lenght-of-stay calculated based on github data (distribution_icu_world.csv) provided by
// https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01726-3#Sec13
static const vector<double> CRITICAL_DURATION_CDF =
    {0.02395, 0.08514, 0.15554, 0.23370, 0.31259, 0.39308, 0.46988, 0.54266, 0.60962, 0.67088,
     0.72574, 0.77359, 0.81612, 0.85151, 0.88142, 0.90602, 0.92574, 0.94139, 0.95484, 0.96481,
     0.97280, 0.97898, 0.98391, 0.98735, 0.99029, 0.99231, 0.99374, 0.99496, 0.99584, 0.99670,
     0.99723, 0.99767, 0.99810, 0.99848, 0.99869, 0.99893, 0.99915, 0.99924, 0.99939, 0.99947,
     0.99955, 0.99962, 0.99969, 0.99977, 0.99980, 0.99984, 0.99986, 0.99987, 0.99990, 0.99990,
     0.99993, 0.99995, 0.99995, 0.99995, 0.99995, 0.99995, 0.99996, 0.99997, 0.99998, 1.00000};


static const float NON_ICU_CRITICAL_MORTALITY = 0.9;          // probability of dying when becoming critical if intensive care is not received

// from Person.h
static const int NUM_AGE_CLASSES = 121;                       // maximum age+1 for a person


struct DynamicParameter {
    DynamicParameter(){};
    DynamicParameter(int s, int d, double v) : start(s), duration(d), value(v) {};
    int start;
    int duration;
    double value;
};


class StrainPars {
  public:
    StrainPars() :
        type(WILDTYPE),
        relInfectiousness(1.0),
        relPathogenicity(1.0),
        relSeverity(1.0),
        relCriticality(1.0),
        relMortality(1.0),
        relIcuMortality(1.0),
        immuneEscapeProb(0.0),
        symptomaticInfectiousPeriod(SYMPTOMATIC_INFECTIOUS_PERIOD),
        relSymptomOnset(1.0) {}
    StrainPars(StrainType t) : StrainPars() { type = t; }
    StrainType type;
    double relInfectiousness;
    double relPathogenicity;
    double relSeverity;
    double relCriticality;
    double relMortality;
    double relIcuMortality;
    double immuneEscapeProb;
    int    symptomaticInfectiousPeriod;
    double relSymptomOnset;
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


/*struct CatchupVaccinationEvent {
    CatchupVaccinationEvent(){};
    CatchupVaccinationEvent(size_t cs, size_t cd, size_t a, double c): campaignStart(cs), campaignDuration(cd), age(a), coverage(c) {};
    size_t campaignStart;
    size_t campaignDuration;
    size_t age;
    double coverage;
};*/


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

    static size_t sampler (const std::vector<double> CDF, const double rand, size_t index = 0) {
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
        const double icu_time_to_death_p = 0.2;
        const size_t icu_time_to_death_n = 2;
        return gsl_ran_negative_binomial(RNG, icu_time_to_death_p, icu_time_to_death_n);
    }

    int sampleCommunityTimeToDeath() const {
        const double daily_prob_of_death = 0.5;
        return gsl_ran_geometric(RNG, daily_prob_of_death);
    }

    double timedInterventionEffect(TimedIntervention ti, size_t day) const;

    unsigned long int randomseed;
    size_t runLength;
    double household_transmission_haz_mult;                      // per-day hazard multiplier of transmission between a co-habitating dyad
    double social_transmission_haz_mult;                         // per-day hazard multiplier of transmission, scaled by the fraction of infectious social contacts
    double workplace_transmission_haz_mult;                      // per-day hazard multiplier of transmission, scaled by the fraction of infectious co-workers
    double school_transmission_haz_mult;                         // per-day hazard multiplier of transmission, scaled by the fraction of infectious students/staff
    double hospital_transmission_haz_mult;                       // per-day hazard multiplier of transmission, scaled by the fraction of infectious staff/patients
    double nursinghome_transmission_haz_mult;                    // per-day hazard multiplier of transmission, scaled by the fraction of infectious staff/residents
    std::vector<double> seasonality;                        // transmissibility multiplier, index by simulation day
    double seasonality_on (const Date *date) const;

    vector<float> susceptibilityByAge;                      // probability of infection given exposure, index by year of age
    vector<float> pathogenicityByAge;                       // probability of clinical disease given infection, index by year of age
    vector<float> severeFractionByAge;                      // probability of severe disease given clinical disease, index by year of age
    map<OutcomeType, vector< vector<float>>> probSeriousOutcome;        // look up probability of severe, critical, fatal outcomes by age & comorbid status
                                                            // vaccine efficacies, indexed by dose
    map<StrainType, vector<double>> VES;                    // vaccine efficacy for susceptibility (can be leaky or all-or-none)
    map<StrainType, vector<double>> VES_NAIVE;              // VES for initially immunologically naive people
    map<StrainType, vector<double>> VEP;                    // vaccine efficacy for pathogenicity
    map<StrainType, vector<double>> VEH;                    // vaccine efficacy against hospitalization, given infection
    map<StrainType, vector<double>> VEF;                    // vaccine efficacy against death, given infection
    map<StrainType, vector<double>> VEI;                    // vaccine efficacy to reduce infectiousness

    int vaccine_dose_to_protection_lag;                     // number of days between when vaccine dose is administered and when its protection begins

    //double IES;                                             // represented using individual-based startingNaturalEfficacy and par->remainingEfficacy
    double IEP;                                             // prior infection efficacy for pathogenicity
    double IEH;                                             // prior infection efficacy against hospitalization, given infection
    double IEF;                                             // prior infection efficacy against death, given infection
    double IEI;                                             // prior infection efficacy to reduce infectiousness

    template<typename T>
    inline T stretchy_vector (const vector<T> &data, size_t idx) const {
        idx = idx < 0 ? 0 : min(idx, data.size() - 1);      // shift idx to be within valid range if outside
        return data.at(idx);                                // using at() to force bounds checking, in case data vector is empty
    }

    double VES_at(size_t dose, StrainType strain = WILDTYPE)       const { return stretchy_vector(VES.at(strain), dose);}
    double VES_NAIVE_at(size_t dose, StrainType strain = WILDTYPE) const { return stretchy_vector(VES_NAIVE.at(strain), dose);}
    double VEP_at(size_t dose, StrainType strain = WILDTYPE)       const { return stretchy_vector(VEP.at(strain), dose);}
    double VEH_at(size_t dose, StrainType strain = WILDTYPE)       const { return stretchy_vector(VEH.at(strain), dose);}
    double VEF_at(size_t dose, StrainType strain = WILDTYPE)       const { return stretchy_vector(VEF.at(strain), dose);}
    double VEI_at(size_t dose, StrainType strain = WILDTYPE)       const { return stretchy_vector(VEI.at(strain), dose);}

    int beginContactTracing;                                // what sim day should contact tracing start
    double contactTracingCoverage;                          // fraction of detected infections that are contact traced as primary cases
    vector<double> contactTracingEV;                        // indexed by LocationType; Expected number of recalled contacts; HOME is used for neighbor contacts
    size_t contactTracingDepth;                             // tracing contacts = 1; tracing contacts-of-contacts = 2; etc...
    int urgent_vax_dose_threshold;                       // what doses can be urgently provided (doses less than this threshold); must be <= numVaccineDoses

    size_t symptom_onset(StrainType strain = WILDTYPE) const { // aka incubation period
        // TODO - CABP: this could be a negative binomial instead (which is discrete)
        double deviate = gsl_ran_gamma(RNG, WILDTYPE_SYMPTOM_ONSET_GAMMA_SHAPE, strainPars[strain].relSymptomOnset * WILDTYPE_SYMPTOM_ONSET_GAMMA_SCALE);
        assert(deviate >= 0);
        return 1 + floor(deviate);
    }
    size_t infectiousness_onset(size_t incubation_pd) const {
        assert(incubation_pd >= 1);
        return gsl_ran_binomial(RNG, INFECTIOUSNESS_ONSET_FRACTION, incubation_pd - 1) + 1;
    }
    size_t pre_severe_symptomatic() const {
        return 1 + floor(gsl_ran_gamma(RNG, PRE_SEVERE_SYMPTOMATIC_GAMMA_SHAPE, PRE_SEVERE_SYMPTOMATIC_GAMMA_SCALE));
    }
    size_t symptom_duration_mild() const {
        return 1 + floor(gsl_ran_gamma(RNG, WILDTYPE_SYMPTOM_DURATION_MILD_GAMMA_SHAPE, WILDTYPE_SYMPTOM_DURATION_MILD_GAMMA_SCALE));
    }
    size_t severe_only_duration() const { return sampler(SEVERE_ONLY_DURATION_CDF, gsl_rng_uniform(RNG)); }
    size_t critical_duration()    const { return sampler(CRITICAL_DURATION_CDF, gsl_rng_uniform(RNG)); }

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

    size_t deathReportingLag (const gsl_rng* REPORTING_RNG) const {
        //return meanDeathReportingLag;
        // neg binom fit to Oct/Nov 2020 data, with mean of ~12 days
        const double n = 2.2786932;
        const double p = 0.1578952;
        return gsl_ran_negative_binomial(REPORTING_RNG, p, n);
    }

    double pathogenicityReduction;                          // == percentage of missed infections in empirical pathogenicity estimates
                                                            // used for interpreting empirical pathogenicity data, *not within simulator*
    double susceptibilityCorrection;                        // susceptibilityCorrection of 0 --> published value; 1 --> everyone 100% susceptible

    void createIcuMortalityReductionModel(double maximum_val, double inflection_sim_day, double slope);
    double icuMortality(ComorbidType comorbidity, size_t age, size_t sim_day) const;
    double icuMortalityFraction;                            // fraction of all deaths that occur in ICUs; Pr{ICU|death}, *NOT* Pr{death|ICU}
                                                            // used for interpreting empirical mortality data, *not within simulation*

                                                            // sign of float determined by sizes of initial and final values
    void createDetectionModel(const vector<vector<double>>& vals, const vector<vector<int>>& inflection_sim_day, const vector<vector<double>>& slopes);
    std::vector<std::vector<double>> probFirstDetection;    // Probability of being *first* detected while {asymp, mild, severe, crit, dead}
                                                            // indexed by day and outcome severity


    //size_t meanDeathReportingLag;                           // number of days from when death occurs to when it's reported
    void createSocialDistancingModel(std::string filename, size_t metric_col, float mobility_logit_shift, float mobility_logit_stretch);
    bool immunityLeaky;                                     // if false, natural and vaccine immunity is all-or-none
    bool retroactiveMatureVaccine;                          // if true, infection causes leaky vaccine to jump from naive to mature protection
    double seroTestFalsePos;                                // probability that seroneg person tests positive -- leaky test
    double seroTestFalseNeg;                                // probability that seropos person tests negative -- leaky test
    size_t numInitialExposed;                               // for general use
    size_t numInitialInfected;                              // for R0 estimation
    double probInitialExposure;                             // for general use, to scale intial exposures with pop size
//    double probInitialInfection;
//    std::vector<double> numDailyExposed;                    // dimension is days
    std::vector<double> probDailyExposure;                  // per person, per day probability of exposure

    std::string populationFilename;
    std::string comorbidityFilename;
    std::string locationFilename;
    std::string publicActivityFilename;
    std::string networkFilename;
    std::string peopleOutputFilename;
    std::string yearlyPeopleOutputFilename;
    std::string dailyOutputFilename;
    std::string rCaseDeathFilename;
    std::string vaccination_file;
    std::string dose_file;
    //std::string annualIntroductionsFilename;                // time series of some external factor determining introduction rate
    //std::vector<double> annualIntroductions;
    //double annualIntroductionsCoef;                         // multiplier to rescale external introductions to something sensible

    double sampleStartingNaturalEfficacy(const gsl_rng* RNG) const {
        // neutralization level relative to mean natural immunity following infection, after: https://www.nature.com/articles/s41591-021-01377-8/figures/1
        const double neut_lvl = pow(2.0, gsl_ran_gaussian(RNG, 1.0));
        const double AESP = covid::util::logistic(3.097703*(log10(neut_lvl/0.2010885)));
        const double AEP = 0.75;
        //const double VEP = 1 - (1 - 0.93)/(1 - 0.52); // values from https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)00675-9/fulltext
        return 1 - (1 - AESP)/(1 - AEP);
    }

    const double __IMMUNE_DECAY_SLOPE   = 0.2;//0.4148295;
    const double __IMMUNE_DECAY_INTCPT  = 2.4603816;

    double _immunityEffectiveStartingTime(double starting_efficacy) const {
        return (log(1.0/starting_efficacy - 1.0) + __IMMUNE_DECAY_INTCPT) / __IMMUNE_DECAY_SLOPE;     // this is a logit transformation
    }

    double remainingEfficacy(double starting_efficacy, double time_delta) const {
        if (not immunityWanes) { return starting_efficacy; }
        const double month_delta = time_delta / 30.0; // published data worked in terms of months
        // KBT fitted logistic model to Pfizer waning data from: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02183-8/fulltext
        //logistic_wane = function(x) {1 / (1 + exp(-(2.4603816 - 0.4148295 * x)))}
        //logit_wane = function(x) { (log(1/x - 1) + 2.4603816)/0.4148295 }
        // first determine where this person's immunity started along the logistic curve
        const double effective_starting_time = _immunityEffectiveStartingTime(starting_efficacy);
        const double l = __IMMUNE_DECAY_INTCPT - __IMMUNE_DECAY_SLOPE * (effective_starting_time + month_delta);
        return covid::util::logistic(l);
    }

    int immunityDuration(const double quantile, const double starting_efficacy) const {
        // logit to determine offset
        const double offset = _immunityEffectiveStartingTime(0.5) - _immunityEffectiveStartingTime(starting_efficacy);
        // sample from logistic distr and add offset to find time to immunity loss
        double immunityDuration = gsl_cdf_logistic_Pinv(quantile, -1.0/__IMMUNE_DECAY_SLOPE) + offset;     // duration in months
        if (immunityDuration < 0) { immunityDuration = 0; }
        //cerr << starting_efficacy << ' ' << offset << ' ' << immunityDuration << endl;
        return 30 * immunityDuration;
    }

    int immunityDuration(const double quantile, const size_t dose, StrainType strain) const {
        const double starting_efficacy = VES_at(dose, strain);
        return immunityDuration(quantile, starting_efficacy);
    }
//    size_t sampleDaysImmune(const gsl_rng* RNG) const {
//        const double t_half = 103;             // days
//        const double threshold = 215;          // below this threshold, susceptible again.  CABP estimated using PHE data and Science paper below
//        const double antibody_init_mean = 3.0; // Rough estimate from https://science.sciencemag.org/content/early/2021/01/06/science.abf4063/
//        const double antibody_init_sd = 0.5;   // ''
//        const double time_to_susceptible = -t_half * log(threshold / pow(10, gsl_ran_gaussian(RNG, antibody_init_sd) + antibody_init_mean)) / log(2);
//        return time_to_susceptible < 0 ? 0 : (size_t) round(time_to_susceptible);
//    }
    bool immunityWanes;                                     // innate waning of both natural and vaccine immunity
    double seroPositivityThreshold;                         // threshold at which an individual would be considered seropositive
    int vaccineImmunityDuration;
    bool vaccineBoosting;                                   // Are we re-vaccinated, either because of waning or because of multi-dose vaccine
    int numVaccineDoses;                                    // Number of times to boost; default is INT_MAX
    std::vector<int> vaccineDoseInterval;                   // How often to we re-vaccinate for initial vaccine course, in days
    int vaccineBoostingInterval;                            // How often to we re-vaccinate for boosting, in days
//    std::vector<CatchupVaccinationEvent> catchupVaccinationEvents;
    int vaccineTargetAge;
    double vaccineTargetCoverage;
    int vaccineTargetStartDate;

    std::vector<StrainPars> strainPars;
    std::vector<std::vector<bool>> crossProtectionMatrix;

    size_t numSurveilledPeople;

    // bool vacCampaign_prioritize_first_doses;
    // bool vacCampaign_flexible_queue_allocation;
    // ReactiveVaccinationStrategyType vacCampaign_reactive_strategy;

    bool traceContacts;                                     // should we keep track of who infected whom

                                                            // e.g. school closures, non-essential business closures, social distancing
    std::map<TimedIntervention, std::vector<double>> timedInterventions;
    std::vector<double> icuMortalityReduction;              // time-varying reduction (from 0) in ICU mortality due to improved Tx

    size_t startDayOfYear;
    size_t startJulianYear;
    bool dailyOutput;
    bool periodicOutput;
    int periodicOutputInterval;
    bool weeklyOutput;
    bool monthlyOutput;
    bool yearlyOutput;
    bool abcVerbose;
    unsigned long int serial;

    VaccineSeroConstraint vaccineSeroConstraint;
    CsmhScenario csmhScenario;

    bool behavioral_autotuning;
    size_t tuning_window;
    size_t num_preview_windows;
    std::string autotuning_dataset;

    std::vector<double> quarantineProbability;
    size_t selfQuarantineDuration;

    bool dump_simulation_data;
};

#endif
