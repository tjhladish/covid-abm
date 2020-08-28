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
    reportedFraction = {0.0, 0.2, 0.8, 0.8, 0.8};       // fraction of asymptomatic, mild, severe, critical, and deaths reported
//    numDailyExposed.push_back(0.0);                     // default: no introductions
    probDailyExposure.push_back(0.0);                   // default: no introductions

    symptomToTestLag = 2;
    defaultReportingLag = 10;
    rlm = nullptr;                                       // reporting lag model
    deathReportingLag = 4;
    numInitialExposed  = 0;
    numInitialInfected = 0;

    pathogenicityModel = ORIGINAL_LOGISTIC;

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

    define_susceptibility_and_pathogenicity();
    mmodsScenario = NUM_OF_MMODS_SCENARIOS; // default to no MMODS scenario
}


void Parameters::define_susceptibility_and_pathogenicity() {
    // values from EDT1 of https://www.medrxiv.org/content/10.1101/2020.03.24.20043018v2.full.pdf
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
    bin_upper = {4, 17, 49, 64, 74, 84, NUM_AGE_CLASSES-1};
    vector<float> severe_fraction = {0.003, 0.001, 0.025, 0.074, 0.122, 0.158, 0.172};

    for (size_t i = 0; i < bin_upper.size(); ++i) {
        const size_t upper_age = bin_upper[i];
        severeFractionByAge.resize(upper_age+1, severe_fraction[i]);
    }

//    for (size_t i = 0; i < susceptibilityByAge.size(); ++i) {
//        cerr << i << "\t" << susceptibilityByAge[i] << "\t" << pathogenicityByAge[i] << endl;
//    }
}


double Parameters::icuMortality(size_t sim_day) const {
    const size_t imr_size = icuMortalityReduction.size();
    //assert(imr_size == 0 or sim_day < (signed) imr_size);
    if (imr_size == 0) {
        return ICU_CRITICAL_MORTALITY;
    } else {
        // if it's a day off the end of the vector, use the last value
        double const imr = icuMortalityReduction.size() > sim_day ? icuMortalityReduction[sim_day] : icuMortalityReduction.back();
        return ICU_CRITICAL_MORTALITY * (1.0 - imr);
    }
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

