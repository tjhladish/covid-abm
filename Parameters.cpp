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

void Parameters::define_defaults() {
    serial = 0;
    randomseed = 5489;
    runLength = 100;
    household_transmissibility = 0.15;
    workplace_transmissibility = 0.15;
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
    reportedFraction = {0.0, 0.01, 0.5, 0.8, 1.0};      // fraction of asymptomatic, mild, severe, critical, and deaths reported
    numDailyExposed.push_back(0.0);                     // default: no introductions

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
}


void Parameters::define_susceptibility_and_pathogenicity() {
    // values from EDT1 of https://www.medrxiv.org/content/10.1101/2020.03.24.20043018v2.full.pdf
    vector<size_t> bin_upper = {9, 19, 29, 39, 49, 59, 69, NUM_AGE_CLASSES-1};
    vector<float> susceptibilities = {0.33, 0.37, 0.69, 0.81, 0.74, 0.8, 0.89, 0.77};
    vector<float> pathogenicities = {0.4, 0.25, 0.37, 0.42, 0.51, 0.59, 0.72, 0.76};

    for (size_t i = 0; i < bin_upper.size(); ++i) {
        const size_t upper_age = bin_upper[i];
        susceptibilityByAge.resize(upper_age+1, susceptibilities[i]);
        pathogenicityByAge.resize(upper_age+1, pathogenicities[i]);
    }

//    for (size_t i = 0; i < susceptibilityByAge.size(); ++i) {
//        cerr << i << "\t" << susceptibilityByAge[i] << "\t" << pathogenicityByAge[i] << endl;
//    }
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

