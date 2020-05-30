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
    household_transmissibility = 0.05;
    workplace_transmissibility = 0.05;
    infectiousOnset = 3;                                // days into infection when infectiousness starts
    infectiousDuration = 10;                            // number of days infectious
    VES = 0.7;
    VES_NAIVE = 0.0;
    VEI = 0.0;
    VEP = 0.0;
    VEH = 0.0;
    hospitalizedFraction = {0.0, 0.15, 0.9};            // rough estimates in mex
    vaccineLeaky = false;
    //secondaryTransmission = true;
    populationFilename = "population.txt";
    locationFilename = "locations.txt";
    peopleOutputFilename = "";
    yearlyPeopleOutputFilename = "";
    dailyOutputFilename = "";
    annualIntroductionsFilename = "";                   // time series of some external factor determining introduction rate
    annualIntroductionsCoef = 1;                        // multiplier to rescale external introductions to something sensible
    annualIntroductions = {1.0};
    daysImmune = 365;
    reportedFraction = {0.0, 0.05, 1.0};                // fraction of asymptomatic, mild, and severe cases reported
    numDailyExposed.push_back(0.0);                     // default: no introductions
    //annualSerotypeFilename = "";

    numInitialExposed  = 0;
    numInitialInfected = 0;

    pathogenicityModel = ORIGINAL_LOGISTIC;
    //primaryRelativeRisk= 0.5;

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

