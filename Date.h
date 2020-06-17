#ifndef __DATE_H
#define __DATE_H

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


class Date {
  public:
    Date():_offset(0),_simulation_day(0) {};
    Date(const Parameters* par):_offset(par->startDayOfYear-1),_simulation_day(0) {};

    int offset()             const { return _offset; }
    inline int day()         const { return _simulation_day; }                // [0, ...]
    int julianDay()          const { return ((day() + offset()) % 365) + 1; } // [1, 365]
    int dayOfMonth()         const { return julianMonth() == 1 ?              // [1, {29,30,31}]
                                            julianDay() :
                                            julianDay() - END_DAY_OF_MONTH[julianMonth()-2]; }
    int nDayPeriod(int n)    const { return (int) day()/n; }                  // [0, ...]
    int week()               const { return (int) day()/7; }                  // [0, ...]
    int julianWeek()         const { return (int) ((julianDay()-1)/7) + 1; }  // [1, 53]
    int month()              const { return _month_ct; }                      // [0, ...]
    int julianMonth()        const {                                          // [1, 12]
        vector<int>::const_iterator it;
        // find first month that hasn't ended (hint: it's this month)
        // julianDay()-1 because this isn't upper_or_equal_bound, which would be convenient
        it = upper_bound(END_DAY_OF_MONTH.begin(), END_DAY_OF_MONTH.end(), julianDay()-1);
        return it - END_DAY_OF_MONTH.begin() + 1; // +1 b/c [1, 12], not [0, 11]
    }
    DayOfWeekType dayOfWeek() const { return (DayOfWeekType) ((day() + offset()) % NUM_OF_DAYS_IN_WEEK); } // [0,6] --> [MON,SUN]
    string monthName()       const { return MONTH_NAMES[julianMonth()-1]; }
    int year()               const { return (int) (day()/365); }

    bool endOfPeriod(int n)  const { return (day()+1) % n == 0; }
    bool endOfWeek()         const { return (day()+1) % 7 == 0; }
    bool endOfMonth()        const {
        vector<int>::const_iterator it;
        // find out if today corresponds to a month-end
        it = find(END_DAY_OF_MONTH.begin(), END_DAY_OF_MONTH.end(), julianDay());
        return it != END_DAY_OF_MONTH.end();
    }
    bool startOfYear()       const { return day() % 365 == 0; }     // is it beginning of 365 day period
    bool endOfYear()         const { return (day()+1) % 365 == 0; } // is it end of 365 day period
    bool startOfJulianYear() const { return julianDay() == 1; }     // is it Jan 1
    bool endOfJulianYear()   const { return julianDay() == 365; }   // is it Dec 31

    void setJulianDay (int julian_day) { _offset = julian_day - day() - 1; }

    void increment() {
        if(endOfMonth()) _month_ct++;
        _simulation_day++;
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

  private:
    int _offset;
    int _simulation_day;
    int _month_ct;
};

#endif
