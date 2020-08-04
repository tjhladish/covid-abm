#ifndef __DATE_H
#define __DATE_H
#include "Parameters.h"

enum DayOfWeekType {
    MONDAY,
    TUESDAY,
    WEDNESDAY,
    THURSDAY,
    FRIDAY,
    SATURDAY,
    SUNDAY,
    NUM_OF_DAYS_IN_WEEK
};

static const std::vector<std::string> DAY_NAMES = {"MON", "TUE", "WED", "THU", "FRI", "SAT", "SUN"};
static const std::vector<std::string> MONTH_NAMES = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
static const std::vector<size_t> COMMON_DAYS_IN_MONTH = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
static const std::vector<size_t> COMMON_END_DAY_OF_MONTH = {31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
static const std::vector<size_t> LEAP_DAYS_IN_MONTH = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
static const std::vector<size_t> LEAP_END_DAY_OF_MONTH = {31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366};


class Date {
  public:
    Date():_simulation_day(0),_month_ct(0),_year_ct(0),_julian_day(1),_julian_year(1) {};
    Date(const Parameters* par):_simulation_day(0),_month_ct(0),_year_ct(0),_julian_day(par->startDayOfYear),_julian_year(par->julianYear) {};

    inline int day()         const { return _simulation_day; }                                     // [0, ...]
    size_t julianDay()       const { return _julian_day; }                                         // [1, {365, 366}]
    size_t dayOfMonth()      const { return julianMonth() == 1 ?                                   // [1, {29,30,31}]
                                            julianDay() :
                                            julianDay() - end_day_of_month()[julianMonth()-2]; }
    size_t nDayPeriod(int n)    const { return (int) day()/n; }                                    // [0, ...]
    size_t week()               const { return (int) day()/NUM_OF_DAYS_IN_WEEK; }                  // [0, ...]
    size_t julianWeek()      const { return (int) ((julianDay()-1)/NUM_OF_DAYS_IN_WEEK) + 1; }     // [1, 53]
    size_t month()              const { return _month_ct; }                                        // [0, ...]
    size_t julianMonth()        const {                                                            // [1, 12]
        vector<size_t>::const_iterator it;
        // find first month that hasn't ended (hint: it's this month)
        // julianDay()-1 because this isn't upper_or_equal_bound, which would be convenient
        const vector<size_t> end_dom = end_day_of_month();
        it = upper_bound(end_dom.begin(), end_dom.end(), julianDay()-1);
        return it - end_dom.begin() + 1; // +1 b/c [1, 12], not [0, 11]
    }

    static DayOfWeekType dayOfWeek(int y, int m, int d) {
        // John Conway's "Doomsday Algorithm"
        // static int t[] = {0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4}; // week starting Sunday
        static int t[] = {6, 2, 1, 4, 6, 2, 4, 0, 3, 5, 1, 3};   // week starting Monday
        y -= m < 3;
        return (DayOfWeekType) ((y + y/4 - y/100 + y/400 + t[m-1] + d) % 7);
    }

    DayOfWeekType dayOfWeek() const { return dayOfWeek(julianYear(), julianMonth(), dayOfMonth()); }


    string dayOfWeekName(int y, int m, int d) const { return DAY_NAMES[dayOfWeek(y, m, d)]; }
    string dayOfWeekName()   const { return DAY_NAMES[dayOfWeek()]; }
    string monthName()       const { return MONTH_NAMES[julianMonth()-1]; }
    size_t year()            const { return _year_ct; }

    bool endOfPeriod(int n)  const { return (day()+1) % n == 0; }
    bool endOfWeek()         const { return (day()+1) % 7 == 0; }
    bool endOfMonth()        const {
        vector<size_t>::const_iterator it;
        // find out if today corresponds to a month-end
        const vector<size_t> end_dom = end_day_of_month();
        it = find(end_dom.begin(), end_dom.end(), julianDay());
        return it != end_dom.end();
    }
    vector<size_t> end_day_of_month() const { return isLeap() ? LEAP_END_DAY_OF_MONTH : COMMON_END_DAY_OF_MONTH; }
    bool startOfYear()       const { return startOfJulianYear(); }
    bool endOfYear()         const { return endOfJulianYear(); }                               // is it end of {365,366} day period
    bool startOfJulianYear() const { return julianDay() == 1; }                                // is it Jan 1
    bool endOfJulianYear()   const { return julianDay() == num_days_in_year(julianYear()); }   // is it Dec 31

    void setJulianDay (int jd) { _julian_day = jd; }
    void setJulianYear (int julian_year) { _julian_year = julian_year; }

    size_t julianYear () const { return _julian_year; }

    void increment() {
        if(endOfMonth()) _month_ct++;
        _simulation_day++;
        _julian_day++;
        if (_julian_day > num_days_in_year(_julian_year)) {
            _julian_day = 1;
            _julian_year++;
        }
    }

    void print() {
        cerr << day() << "\t" << julianDay() << "\t" << year()
             << "\t(" << monthName() << " " << dayOfMonth() << ")\t" << month() << "\t" << julianMonth();
        if (endOfWeek()) cerr << " EoW";
        if (endOfMonth()) cerr << " EoM";
        if (startOfYear()) cerr << " SoY";
        if (endOfYear()) cerr << " EoY";
        if (endOfJulianYear()) cerr << " EoJY";
        cerr << endl;
    }

    bool isLeap() const {
        return Date::isLeap(julianYear());
    }

    static bool isLeap (size_t julian_year) {
        if (julian_year % 4 != 0) {
            return false;
        } else if (julian_year % 100 != 0) {
            return true;
        } else if (julian_year % 400 != 0) {
            return false;
        } else {
            return true;
        }
    }

    static size_t num_days_in_year(size_t julian_year) {
        if (isLeap(julian_year)) {
            return LEAP_END_DAY_OF_MONTH.back();
        } else {
            return COMMON_END_DAY_OF_MONTH.back();
        }
    }

    string to_string(vector<string> format, string sep="/", bool zero_pad = true) const {
        stringstream ss;
        for (size_t i = 0; i<format.size(); ++i) {
            string part = format[i];
            if (zero_pad) { ss << setfill('0'); }
            if (part == "yyyy") {
                ss << setw(4) << julianYear();
            } else if (part == "mm") {
                ss << setw(2) << julianMonth();
            } else if (part == "dd") {
                ss << setw(2) << dayOfMonth();
            } else {
                cerr << "ERROR: Date::date_string() format argument \"" << part << "\" not supported.  Please use yyyy, mm, and/or dd.\n";
                exit (-1);
            }
            if (i < format.size() - 1) { ss << sep; }
        }
        return ss.str();
    }

  private:
    int _simulation_day;
    size_t _month_ct;
    size_t _year_ct;
    size_t _julian_day;
    size_t _julian_year;
};

#endif
