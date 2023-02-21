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
// #include "Vac_Campaign.h"

using namespace covid::util;

void Parameters::define_defaults() {
    serial                         = 0;
    randomseed                     = 5489;
    runLength                      = 0;

    household_transmission_haz_mult   = 0.15;
    social_transmission_haz_mult      = 0.15;
    workplace_transmission_haz_mult   = 0.15;
    school_transmission_haz_mult      = 0.15;
    hospital_transmission_haz_mult    = 0.015;
    nursinghome_transmission_haz_mult = 0.3;

    VES                            = {{WILDTYPE, {0.7}}};
    VES_NAIVE                      = {{WILDTYPE, {0.0}}};
    VEP                            = {{WILDTYPE, {0.0}}};
    VEH                            = {{WILDTYPE, {0.0}}};
    VEF                            = {{WILDTYPE, {0.0}}};
    VEI                            = {{WILDTYPE, {0.0}}};

    //IES                            = 0.0;
    IEP                            = 0.75;
    IEH                            = 0.0;
    IEF                            = 0.0;
    IEI                            = 0.0;

    vaccine_dose_to_protection_lag = 10;
    immunityLeaky                  = false;
    beginContactTracing            = INT_MAX;
    contactTracingCoverage         = 0.0;
    contactTracingEV               = vector<double>(NUM_OF_LOCATION_TYPES);
    contactTracingDepth            = 0;
    quarantineProbability          = vector<double>(contactTracingDepth);
    quarantineDuration             = 0;
    urgent_vax_dose_threshold      = 0;  // default to 0 (no urgent vaccines) but must be <= numVaccineDoses
    populationFilename             = "population.txt";
    locationFilename               = "locations.txt";
    publicActivityFilename         = "public-activity.txt";
    networkFilename                = "network.txt";
    vaccinationFilename            = "";
    doseFilename                   = "";
    riskGroupsFilename             = "";

    behaviorInputFilename               = "autotuned-behavior.csv";
    behaviorOutputFilename             = "autotuned-behavior.csv";
    rCaseDeathFilename             = "";

    peopleOutputFilename           = "";
    yearlyPeopleOutputFilename     = "";
    dailyOutputFilename            = "";
    probDailyExposure.push_back(0.0);                   // default: no introductions
    icuMortalityFraction           = 0.5;               // fraction of empirical deaths that are assumed to have occured in ICUs
    pathogenicityReduction         = 0.0;               // fraction of empirical infections that were missed in pathogenicity studies
    susceptibilityCorrection       = 0.0;               // 0.0 means use published age-specific susceptibility values; 1.0 means 100% susceptibility

    symptomToTestLag               = 2;
    defaultReportingLag            = 10;
    rlm                            = nullptr;           // reporting lag model
    numInitialExposed              = 0;
    numInitialInfected             = 0;
    probInitialExposure            = 0.0;
    //probInitialInfection           = 0.0;

    numVaccineDoses                = 3;
    vaccineDoseInterval            = vector<int>(numVaccineDoses);

    immunityWanes                  = false;
    seropositivityThreshold        = 0.0;
    seroconversionLag              = 10.0;
    vaccineImmunityDuration        = INT_MAX;
    vaccineBoosting                = false;
    vaccineBoostingInterval        = 730;
    retroactiveMatureVaccine       = false;

    // vacCampaign_prioritize_first_doses = false;
    // vacCampaign_flexible_queue_allocation = false;
    // vacCampaign_reactive_strategy = NUM_OF_REACTIVE_VAC_STRATEGY_TYPES;

    for (int strain = 0; strain < NUM_OF_STRAIN_TYPES; ++strain) {
        strainPars.emplace_back((StrainType) strain);
    }

    numSurveilledPeople            = INT_MAX;

    traceContacts                  = false;
    startDayOfYear                 = 1;

    dailyOutput                    = false;
    periodicOutput                 = false;
    periodicOutputInterval         = 5;
    weeklyOutput                   = false;
    monthlyOutput                  = false;
    yearlyOutput                   = false;
    abcVerbose                     = false;

    vaccineInfConstraint           = VACCINATE_ALL_INF_STATUSES;
    vaccineSeroConstraint          = VACCINATE_ALL_SERO_STATUSES;
    seroTestFalsePos               = 0.0;
    seroTestFalseNeg               = 0.0;

//    csmhScenario = NUM_OF_CSMH_SCENARIOS; // default to no scenario

    behavioral_autotuning          = false;
    behavior_fitting_data_target   = NUM_OF_AUTO_FITTING_DATA_TARGETS;
    death_tuning_offset            = 0;
    tuning_window                  = INT_MAX;
    num_preview_windows            = INT_MAX;
    // autotuning_dataset             = "";
    dump_simulation_data           = false;
}


void Parameters::define_susceptibility_and_pathogenicity() {
    // values from Extended Data Fig. 4 of
    // https://www.nature.com/articles/s41591-020-0962-9#Sec12
    // now published in Nat Med
    vector<size_t> bin_upper = {9, 19, 29, 39, 49, 59, 69, NUM_AGE_CLASSES-1};
//    vector<float> susceptibilities(8, 1.0); // made up values
    //                                   9,   19,   29,   39,   49,  59,   69,   120
   // vector<float> susceptibilities = {0.40, 0.38, 0.79, 0.86, 0.80, 0.82, 0.88, 0.74};
    vector<float> susceptibilities = vector<float>(8, 1.0);
    vector<float> pathogenicities  = {0.29, 0.21, 0.27, 0.33, 0.40, 0.49, 0.63, 0.69};

    for (size_t i = 0; i < bin_upper.size(); ++i) {
        const size_t upper_age = bin_upper[i];
        // susceptibilityCorrection of 0 --> published value; 1 --> 100% susceptible
        //susceptibilityByAge.resize(upper_age+1, 1.0 - (1.0 - susceptibilityCorrection)*(1.0 - susceptibilities[i]));
        //pathogenicityByAge.resize(upper_age+1, (1.0 - pathogenicityReduction)*pathogenicities[i]);
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

    // icu admissions (out of entire cohort, not just those hospitalized) -- this comment seems contradictory with comment on 153
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
//cerr << "bin max age, healthy, comorbid icu:hosp ratio: " << bin_upper[i] << " " << critical_com_neg[i]/severe_com_neg[i] << " " << critical_com_pos[i]/severe_com_pos[i] << endl;
    }

    // numbers for people over age 59 suggest a majority die outside of ICU
    // the model separately handles deaths outside of ICU (with higher mortality)
    bin_upper = {9, 19, 29, 39, 49, NUM_AGE_CLASSES-1};
    for (size_t i = 0; i < bin_upper.size(); ++i) {
        const size_t upper_age = bin_upper[i];
        probSeriousOutcome[DEATH][HEALTHY ].resize(upper_age+1, icuMortalityFraction*fatal_com_neg[i]/critical_com_neg[i]);
        probSeriousOutcome[DEATH][COMORBID].resize(upper_age+1, icuMortalityFraction*fatal_com_pos[i]/critical_com_pos[i]);
    }
}


void Parameters::createIcuMortalityReductionModel(double maximum_val, double inflection_sim_day, double slope) {
    icuMortalityReduction = vector<double>(runLength);
    for (size_t sim_day = 0; sim_day < runLength; ++sim_day) {
        icuMortalityReduction[sim_day] = maximum_val * logistic( slope*((double) sim_day - inflection_sim_day) );
        //cerr << "mortality reduction: " << Date::to_ymd(2020, sim_day + startDayOfYear) << " " << icuMortalityReduction[sim_day] << endl;
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


void Parameters::createDetectionModel(const vector<vector<double>>& vals, const vector<vector<int>>& inflection_sim_day, const vector<vector<double>>& slopes) {
    assert(inflection_sim_day.size() == (vals.size() - 1) and slopes.size() == (vals.size() - 1));
    for (const vector<double>& v : vals)            { assert(v.size() == NUM_OF_OUTCOME_TYPES); }
    for (const vector<int>& v : inflection_sim_day) { assert(v.size() == NUM_OF_OUTCOME_TYPES); }
    for (const vector<double>& v : slopes)          { assert(v.size() == NUM_OF_OUTCOME_TYPES); }

    for (size_t outcome = 0; outcome < NUM_OF_OUTCOME_TYPES; ++outcome) {
        for (const vector<double>& v : vals) {
            if (v[outcome] < 0.0 or v[outcome] > 1.0) { cerr << "WARNING: Detection probability out-of-bounds for outcome type: " << (OutcomeType) outcome << " [" << v[outcome] << "]" << endl; }
        }
    }

    probFirstDetection = vector<vector<double>>(runLength, vector<double>(NUM_OF_OUTCOME_TYPES, 0.0));
    for (size_t sim_day = 0; sim_day < runLength; ++sim_day) {
        vector<double> probs(NUM_OF_OUTCOME_TYPES);
        for (size_t outcome = 0; outcome < NUM_OF_OUTCOME_TYPES; ++outcome) {
            probs[outcome] = vals[0][outcome];
            for (size_t i = 1; i < vals.size(); ++i) {
                probs[outcome] += (vals[i][outcome] - vals[i-1][outcome]) * logistic( slopes[i-1][outcome] * ((double) sim_day - inflection_sim_day[i-1][outcome]) );
            }
//if (outcome == 0) cerr << "d, [i,m,f], val: " << (int) sim_day << " [" << vals[0][outcome] << ", " << vals[1][outcome] << ", " << vals[2][outcome] << "] " << probs[outcome] << endl;
        }
        probFirstDetection[sim_day] = probs;
    }
// These next five lines will output the detection probabilities over time
//    cerr << "day asymp mild severe crit death\n";
//    for (unsigned int i = 0; probFirstDetection.size(); ++i) {
//        cerr << i << " "; cerr_vector(toReportedFraction(probFirstDetection[i])); cerr << endl;
//    }
//    exit(1);
}


void Parameters::createReportingLagModel(std::string filename) { rlm = new ReportingLagModel(filename); }

double Parameters::seasonality_on (const Date *date) const { return seasonality.at(date->julianDay() - 1); }

void Parameters::createSocialDistancingModel(std::string filename, size_t metric_col, float mobility_logit_shift, float mobility_logit_stretch) {
    // expects that the mobility data being read in has a first column with increasing, consecutive dates
    // and a second column with the mobility numbers (on [0,1])
    timedInterventions[SOCIAL_DISTANCING].clear();
    //timedInterventions[SOCIAL_DISTANCING].resize(run_length, 0.0);

    const string sim_start_date = Date::to_ymd(startJulianYear, startDayOfYear);

    std::vector<std::vector<std::string>> mobility_data = covid::util::read_2D_vector_file(filename, ',');
    bool header = true;
    for (size_t i = header; i < mobility_data.size(); ++i) {
        vector<string> row = mobility_data[i];
        assert(row.size() >= 2);
        const string date     = row[0];
        const double mobility = stod(row[metric_col]);
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


size_t ReportingLagModel::sample(const gsl_rng* REPORTING_RNG, const Date* date) const {
    return sample(REPORTING_RNG, date->to_string({"yyyy", "mm", "dd"}, "-"));
}


double Parameters::timedInterventionEffect(TimedIntervention ti, size_t day) const {
    return timedInterventions.at(ti)[day];
}


//void Parameters::loadAnnualIntroductions(string annualIntrosFilename) {
//    ifstream iss(annualIntrosFilename.c_str());
//    if (!iss) {
//        cerr << "ERROR: " << annualIntrosFilename << " not found." << endl;
//        exit(114);
//    }
//    annualIntroductions.clear();
//
//    char buffer[500];
//    double intros;
//    istringstream line(buffer);
//
//    while (iss) {
//        iss.getline(buffer,500);
//        line.clear();
//        line.str(buffer);
//        if (line >> intros) {
//            annualIntroductions.push_back(intros);
//        }
//    }
//
//    iss.close();
//    return;
//}
