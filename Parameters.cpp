#include "Parameters.h"
#include "Location.h"
#include "Utility.h"
#include <fstream>
#include <sstream>
#include <numeric>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <climits> // INT_MAX
#include "Date.h"

using namespace covid::util;

void Parameters::define_defaults() {
    serial = 0;
    randomseed = 5489;
    runLength = 0;
    household_transmissibility = 0.15;
    workplace_transmissibility = 0.15;
    school_transmissibility = 0.15;
    social_transmissibility = 0.15;
    VES = 0.7;
    VES_NAIVE = 0.0;
    VEI = 0.0;
    VEP = 0.0;
    VEH = 0.0;
    vaccineLeaky = false;
    //secondaryTransmission = true;
    populationFilename = "population.txt";
    locationFilename = "locations.txt";
    networkFilename = "network.txt";
    peopleOutputFilename = "";
    yearlyPeopleOutputFilename = "";
    dailyOutputFilename = "";
    annualIntroductionsFilename = "";                   // time series of some external factor determining introduction rate
    annualIntroductionsCoef = 1;                        // multiplier to rescale external introductions to something sensible
    annualIntroductions = {1.0};
    daysImmune = 365;
    probFirstDetection = {0.0, 0.1, 0.6, 0.3, 0.1};       // prob of detection if not detected earlier {asymptomatic, mild, severe, critical, deaths}
//    numDailyExposed.push_back(0.0);                     // default: no introductions
    probDailyExposure.push_back(0.0);                   // default: no introductions
    icuMortalityFraction = 0.5;                         // fraction of empirical deaths that are assumed to have occured in ICUs

    symptomToTestLag = 2;
    defaultReportingLag = 10;
    rlm = nullptr;                                       // reporting lag model
    deathReportingLag = 4;
    numInitialExposed  = 0;
    numInitialInfected = 0;
    probInitialExposure = 0.0;
    //probInitialInfection = 0.0;

    catchupVaccinationEvents.clear();
    vaccineTargetAge = 9;
    vaccineTargetCoverage = 0.0;
    vaccineTargetStartDate = INT_MAX;
    numVaccineDoses = 3;
    vaccineDoseInterval = 182;

    linearlyWaningVaccine = false;
    vaccineImmunityDuration = INT_MAX;
    vaccineBoosting = false;
    vaccineBoostingInterval = 730;
    retroactiveMatureVaccine = false;

    numSurveilledPeople = INT_MAX;

    traceContacts = false;
    startDayOfYear = 1;

    dailyOutput   = false;
    periodicOutput  = false;
    periodicOutputInterval  = 5;
    weeklyOutput  = false;
    monthlyOutput = false;
    yearlyOutput  = false;
    abcVerbose    = false;

    // WHO vaccine mechanism variables
    vaccineSeroConstraint = VACCINATE_ALL_SERO_STATUSES;
    seroTestFalsePos = 0.0;
    seroTestFalseNeg = 0.0;

    mmodsScenario = NUM_OF_MMODS_SCENARIOS; // default to no MMODS scenario
}


void Parameters::define_susceptibility_and_pathogenicity() {
    // values from Extended Data Fig. 4 of
    // https://www.nature.com/articles/s41591-020-0962-9#Sec12
    // now published in Nat Med
    vector<size_t> bin_upper = {9, 19, 29, 39, 49, 59, 69, NUM_AGE_CLASSES-1};
    vector<float> susceptibilities = {0.33, 0.37, 0.69, 0.81, 0.74, 0.8, 0.89, 0.77};
    vector<float> pathogenicities = {0.4, 0.25, 0.37, 0.42, 0.51, 0.59, 0.72, 0.76};

    for (size_t i = 0; i < bin_upper.size(); ++i) {
        const size_t upper_age = bin_upper[i];
        susceptibilityByAge.resize(upper_age+1, susceptibilities[i]);
        pathogenicityByAge.resize(upper_age+1, pathogenicities[i]);
    }

    // https://www.cdc.gov/mmwr/volumes/69/wr/mm6915e3.htm
    //bin_upper = {4, 17, 49, 64, 74, 84, NUM_AGE_CLASSES-1};
    //vector<float> severe_fraction = {0.003, 0.001, 0.025, 0.074, 0.122, 0.158, 0.172};

    //for (size_t i = 0; i < bin_upper.size(); ++i) {
    //    const size_t upper_age = bin_upper[i];
    //    severeFractionByAge.resize(upper_age+1, severe_fraction[i]);
    //}

    // https://www.cdc.gov/mmwr/volumes/69/wr/mm6924e2.htm
    // SEVERE, CRITICAL, DEATH
    bin_upper = {9, 19, 29, 39, 49, 59, 69, 79, NUM_AGE_CLASSES-1};

    // based on hospital admissions (including ICU)
    //                                      9,      19,      29,      39,      49,      59,      69,      79,     120
    vector<float> severe_com_neg   = {0.03689, 0.02279, 0.02688, 0.04445, 0.06441, 0.09570, 0.15355, 0.27867, 0.30095};
    vector<float> severe_com_pos   = {0.22294, 0.14884, 0.17505, 0.24209, 0.29597, 0.36328, 0.49908, 0.64716, 0.62338};

    // icu admissions (out of entire cohort, not just those hospitalized)
    vector<float> critical_com_neg = {0.00703, 0.00337, 0.00302, 0.00725, 0.01267, 0.02053, 0.03675, 0.07110, 0.05189};
    vector<float> critical_com_pos = {0.05008, 0.03468, 0.03369, 0.05298, 0.06374, 0.08276, 0.10870, 0.11933, 0.07465};

    vector<float> fatal_com_neg    = {0.00088, 0.00079, 0.00130, 0.00113, 0.00353, 0.00908, 0.02361, 0.10218, 0.29805};
    vector<float> fatal_com_pos    = {0.00646, 0.00771, 0.01370, 0.02767, 0.04458, 0.07837, 0.16704, 0.31670, 0.49668};

    probSeriousOutcome[SEVERE]   = vector<vector<float>>(NUM_OF_COMORBID_TYPES, vector<float>());
    probSeriousOutcome[CRITICAL] = vector<vector<float>>(NUM_OF_COMORBID_TYPES, vector<float>());
    probSeriousOutcome[DEATH]    = vector<vector<float>>(NUM_OF_COMORBID_TYPES, vector<float>());

    for (size_t i = 0; i < bin_upper.size(); ++i) {
        const size_t upper_age = bin_upper[i];
        //probSeriousOutcome[SEVERE  ][HEALTHY ].resize(upper_age+1, 0.5*severe_com_neg[i]);
        //probSeriousOutcome[SEVERE  ][COMORBID].resize(upper_age+1, 0.5*severe_com_pos[i]);
        probSeriousOutcome[SEVERE  ][HEALTHY ].resize(upper_age+1, severe_com_neg[i]);
        probSeriousOutcome[SEVERE  ][COMORBID].resize(upper_age+1, severe_com_pos[i]);
        probSeriousOutcome[CRITICAL][HEALTHY ].resize(upper_age+1, critical_com_neg[i]/severe_com_neg[i]);
        probSeriousOutcome[CRITICAL][COMORBID].resize(upper_age+1, critical_com_pos[i]/severe_com_pos[i]);
    }

    // numbers for people over age 59 suggest a majority die outside of ICU
    // the model separately handles deaths outside of ICU (with mortality = 1.0)
    bin_upper = {9, 19, 29, 39, 49, NUM_AGE_CLASSES-1};
    for (size_t i = 0; i < bin_upper.size(); ++i) {
        const size_t upper_age = bin_upper[i];
        probSeriousOutcome[DEATH][HEALTHY ].resize(upper_age+1, icuMortalityFraction*fatal_com_neg[i]/critical_com_neg[i]);
        probSeriousOutcome[DEATH][COMORBID].resize(upper_age+1, icuMortalityFraction*fatal_com_pos[i]/critical_com_pos[i]);
    }
}


double Parameters::icuMortality(ComorbidType comorbidity, size_t age, size_t sim_day) const {
    double imr = 0.0;
    if (icuMortalityReduction.size() > 0) {
        // if it's a day off the end of the vector, use the last value
        imr = icuMortalityReduction.size() > sim_day ? icuMortalityReduction[sim_day] : icuMortalityReduction.back();
    }
    return probSeriousOutcome.at(DEATH)[comorbidity][age] * (1.0 - imr);
}


void Parameters::createIcuMortalityReductionModel(double maximum_val, double inflection_sim_day, double slope) {
    icuMortalityReduction = vector<double>(runLength);
    for (size_t sim_day = 0; sim_day < runLength; ++sim_day) {
        icuMortalityReduction[sim_day] = maximum_val * logistic( slope*(sim_day - inflection_sim_day) );
        //cerr << "mortality reduction: " << Date::to_ymd(2020, sim_day + startDayOfYear) << " " << icuMortalityReduction[sim_day] << endl;
    }
}


void Parameters::createReportingLagModel(std::string filename) {
    rlm = new ReportingLagModel(filename);
}


void Parameters::createSocialDistancingModel(std::string filename, float mobility_logit_shift, float mobility_logit_stretch) {
    // expects that the mobility data being read in has a first column with increasing, consecutive dates
    // and a second column with the mobility numbers (on [0,1])
    timedInterventions[SOCIAL_DISTANCING].clear();
    //timedInterventions[SOCIAL_DISTANCING].resize(run_length, 0.0);

    const string sim_start_date = Date::to_ymd(julianYear, startDayOfYear);

    std::vector<std::vector<std::string>> mobility_data = covid::util::read_2D_vector_file(filename, ',');
    bool header = true;
    for (size_t i = header; i < mobility_data.size(); ++i) {
        vector<string> row = mobility_data[i];
        assert(row.size() >= 2);
        const string date     = row[0];
        const double mobility = stod(row[1]);
        if (date < sim_start_date) {
            continue;
        } else if (date > sim_start_date and timedInterventions[SOCIAL_DISTANCING].size() == 0) {
            cerr << "ERROR: simulation starts on " << sim_start_date << ", and mobility input file starts too late (on " << date << ")\n";
            exit(19);
        } else {
            const double SD = covid::util::logistic((covid::util::logit(mobility) + mobility_logit_shift) * mobility_logit_stretch);
//cerr << "storing " << date << ", " << mobility << " as " << SD << endl;
            timedInterventions[SOCIAL_DISTANCING].push_back(SD);
        }
    }
    const double last_value = timedInterventions[SOCIAL_DISTANCING].back();
    timedInterventions[SOCIAL_DISTANCING].resize(runLength, last_value);
}


size_t ReportingLagModel::sample(const gsl_rng* RNG, const Date* date) const {
    return sample(RNG, date->to_string({"yyyy", "mm", "dd"}, "-"));
}


double Parameters::timedInterventionEffect(TimedIntervention ti, size_t day) const {
    return timedInterventions.at(ti)[day];
}


void Parameters::loadAnnualIntroductions(string annualIntrosFilename) {
    ifstream iss(annualIntrosFilename.c_str());
    if (!iss) {
        cerr << "ERROR: " << annualIntrosFilename << " not found." << endl;
        exit(114);
    }
    annualIntroductions.clear();

    char buffer[500];
    double intros;
    istringstream line(buffer);

    while (iss) {
        iss.getline(buffer,500);
        line.clear();
        line.str(buffer);
        if (line >> intros) {
            annualIntroductions.push_back(intros);
        }
    }

    iss.close();
    return;
}

