#include <iostream>
#include "Date.h"
#include "Parameters.h"

// compile with e.g.: g++ -O2 -std=c++17 -Wall -Wextra -Wno-deprecated-declarations --pedantic date_test.cpp -o date_test /home/tjhladish/work/covid-abm/Parameters.o /home/tjhladish/work/covid-abm/Utility.o -lgsl

using namespace std;

const gsl_rng* RNG_local = gsl_rng_alloc (gsl_rng_taus2);

int main() {
    for (size_t y = 1998; y < 2010; ++y) {
        //date->setJulianYear(y);
        //cerr << date->julianYear();
        //if (date->isLeap()) { // works!
        cerr << y;
        if (Date::isLeap(y)) {
            cerr << " is a leap year and has " << Date::num_days_in_year(y) << " days\n";
        } else {
            cerr << " is not a leap year and has " << Date::num_days_in_year(y) << " days\n";
        }
    }
    cerr << "Jul 31, 2020: " << Date::dayOfWeek(2020, 7, 31) << endl;
    cerr << "Feb 29, 2020: " << Date::dayOfWeek(2020, 2, 29) << endl;
    cerr << "Feb 20, 2019: " << Date::dayOfWeek(2019, 2, 20) << endl;

    gsl_rng_set(RNG_local, 4);
    ReportingLagModel rlm;
    rlm.read_csv("case_report_delay.csv");
/*    vector<vector<string>> reporting_lag_data = covid::util::read_2D_vector_file("case_report_delay.csv", ',');
    cerr << reporting_lag_data[0][0] << endl;
    bool header = true;
    for (size_t i = header; i < reporting_lag_data.size(); ++i) {
        vector<string> row = reporting_lag_data[i];
        assert(row.size() == 5);
        string date    = row[0];
        double a_shape = stod(row[3]);
        double b_scale = stod(row[4]);
        rlm.insert_lag(date, a_shape, b_scale);
    }*/
    cerr << "earliest lag date: " << rlm.earliest_date() << endl;
    cerr << "latest lag date: " << rlm.latest_date() << endl;
    cerr << "3/1 dev: " << rlm.sample(RNG_local, "2020-02-01") << endl;
    cerr << "4/1 dev: " << rlm.sample(RNG_local, "2020-03-01") << endl;
    cerr << "5/1 dev: " << rlm.sample(RNG_local, "2020-04-01") << endl;
    cerr << "6/1 dev: " << rlm.sample(RNG_local, "2020-05-01") << endl;
    cerr << "7/1 dev: " << rlm.sample(RNG_local, "2020-06-01") << endl;
    cerr << "12/1 dev: " << rlm.sample(RNG_local, "2020-12-01") << endl;

    //static size_t to_julian_day(string date_string) {
    //static string to_ymd (size_t julian_day, int julian_year) {

    cerr << "2020-01-01 jd: " << Date::to_julian_day("2020-01-01") << endl;
    cerr << "2020-02-01 jd: " << Date::to_julian_day("2020-02-01") << endl;
    cerr << "2020-03-01 jd: " << Date::to_julian_day("2020-03-01") << endl;
    cerr << "2020-12-31 jd: " << Date::to_julian_day("2020-12-31") << endl;
    cerr << "2021-12-31 jd: " << Date::to_julian_day("2021-12-31") << endl;

    cerr << "2020, day 1   ymd: "   << Date::to_ymd(2020,   1) << endl;
    cerr << "2020, day 32  ymd: "  << Date::to_ymd(2020,  32) << endl;
    cerr << "2020, day 61  ymd: "  << Date::to_ymd(2020,  61) << endl;
    cerr << "2020, day 366 ymd: " << Date::to_ymd(2020, 366) << endl;
    cerr << "2021, day 365 ymd: " << Date::to_ymd(2021, 365) << endl;

    Parameters* par = new Parameters();
    par->startJulianYear = 2020;
    par->startDayOfYear  = 32;

    cerr << "sim day     0 (2020, julian start = Feb 1): " << Date::to_ymd(    0, par) << endl;
    cerr << "sim day   334 (2020, julian start = Feb 1): " << Date::to_ymd(  334, par) << endl;
    cerr << "sim day   335 (2020, julian start = Feb 1): " << Date::to_ymd(  335, par) << endl;
    cerr << "sim day  1000 (2020, julian start = Feb 1): " << Date::to_ymd( 1000, par) << endl;
    cerr << "sim day   -32 (2020, julian start = Feb 1): " << Date::to_ymd(  -32, par) << endl;
    cerr << "sim day -1000 (2020, julian start = Feb 1): " << Date::to_ymd(-1000, par) << endl;

    Date* date = new Date();
    date->setJulianYear(2020);
    date->setJulianDay(257);
    cerr << "year 2020, day 257: "          << date->to_ymd() << endl;

    cerr << "All of these should be true (1):" << endl;
    cerr << "2020-09-12 is earlier?:\t\t"      << ((string) "2020-09-12" <  *date) << endl;
    cerr << "2020-09-12 is earlier?:\t\t"      << ("2020-09-12" <  *date) << endl;
    cerr << "2020-09-12 is not later?:\t"    << ("2020-09-12" <= *date) << endl;
    cerr << "2020-09-12 is not same day?:\t" << ("2020-09-12" != *date) << endl;
    cerr << "2020-09-13 is same day?:\t"     << ("2020-09-13" == *date) << endl;
    cerr << "2020-09-14 is later?:\t\t"        << ("2020-09-14" >  *date) << endl;
    cerr << "2020-09-14 is not earlier?:\t"  << ("2020-09-14" >= *date) << endl;

    cerr << "2020-09-12 is earlier?:\t\t"      << (*date >  "2020-09-12") << endl;
    cerr << "2020-09-12 is not later?:\t"    << (*date >= "2020-09-12") << endl;
    cerr << "2020-09-12 is not same day?:\t" << (*date != "2020-09-12") << endl;
    cerr << "2020-09-13 is same day?:\t"     << (*date == "2020-09-13") << endl;
    cerr << "2020-09-14 is later?:\t\t"        << (*date <  "2020-09-14") << endl;
    cerr << "2020-09-14 is not earlier?:\t"  << (*date <= "2020-09-14") << endl;

    cerr << "\nincrement and decrement testing" << endl;
    date->setJulianDay(364);
    for (int i = 0; i < 35; ++i) { date->increment(); date->print(); }
    for (int i = 0; i < 35; ++i) { date->print(); date->decrement(); }
}
