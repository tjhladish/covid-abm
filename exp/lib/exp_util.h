#ifndef _EXP_UTIL_H
#define _EXP_UTIL_H

double calculate_conditional_death_reporting_probability(const vector<double> &rho_vals, double RF_death) {
    const double rho_death  = 1.0 - (1.0 - RF_death)/((1.0 - rho_vals[0])*(1.0 - rho_vals[1])*(1.0 - rho_vals[2])*(1.0 - rho_vals[3]));
    return rho_death;
}

void add_death_probabilities(vector<vector<double>>& first_detection_probs, const double RF_death) {
    for (auto && v: first_detection_probs) { v.push_back(calculate_conditional_death_reporting_probability(v, RF_death)); }
}

vector<vector<int>> create_sim_day_matrix(const Parameters* par, vector<string>& inflection_dates) {
    vector<vector<int>> inflection_matrix;
    for (string date_string: inflection_dates) {
        const int sim_day = Date::to_sim_day(par->startJulianYear, par->startDayOfYear, date_string);
        inflection_matrix.push_back(vector<int>(NUM_OF_OUTCOME_TYPES, sim_day)); // assume same inflection date for all outcome types
    }
    return inflection_matrix;
}

vector<vector<double>> create_slope_matrix(size_t num_inflections, double slope) { // assume same logistic slope for all reporting transitions
    vector<vector<double>> slope_matrix;
    for (size_t i = 0; i < num_inflections; ++i) {
        slope_matrix.push_back(vector<double>(NUM_OF_OUTCOME_TYPES, slope));
    }
    return slope_matrix;
}

vector<vector<double>> as_reported_fractions(const Parameters* par, const vector<vector<double>>& first_detection_probs) {
    vector<vector<double>> RF_matrix;
    for (auto&& v: first_detection_probs) {
        RF_matrix.push_back(par->toReportedFraction(v));
    }
    return RF_matrix;
}

void cerr_matrix(const vector<vector<double>> mat) {
    for (auto v: mat) {
        cerr_vector(v); cerr << endl;
    }
}

// Take a list of values, return original indices sorted by value
vector<int> ordered(vector<int> const& values) {

    vector<pair<int,int> > pairs(values.size());
    for(size_t pos=0; pos<values.size(); pos++) {
        pairs[pos] = make_pair(values[pos],pos);
    }

    //bool comparator ( const mypair& l, const mypair& r) { return l.first < r.first; }
    std::sort( pairs.rbegin(), pairs.rend() ); // sort greatest to least
    vector<int> indices(values.size());
    for(size_t i=0; i < pairs.size(); i++) indices[i] = pairs[i].second;

    return indices;
}


string calculate_process_id(vector<double> &args, string &argstring) {
    // CCRC32 checksum based on string version of argument values
    CCRC32 crc32;
    crc32.Initialize();

    for (size_t i = 0; i < args.size(); i++) argstring += to_string((double) args[i]) + " ";

    const unsigned char* argchars = reinterpret_cast<const unsigned char*> (argstring.c_str());
    const int len = argstring.length();
    const int process_id = crc32.FullCRC(argchars, len);

    return to_string(process_id);
}


string report_process_id (vector<double> &args, const unsigned long int serial, const time_t global_start_time, const time_t replicate_start_time) {
    double dif = difftime (replicate_start_time, global_start_time);

    string argstring;
    const string process_id = calculate_process_id(args, argstring);

    cerr << "pid in report_process_id (num args = " << args.size() << "): " << process_id << endl;
    stringstream ss;
    ss << "begin " << process_id << " " << dec << serial << " " << dif << " " << argstring << endl;
    string output = ss.str();
    fputs(output.c_str(), stderr);

    return to_string(process_id);
}


void append_if_finite(vector<double> &vec, double val) {
    if (isfinite(val)) {
        vec.push_back((double) val);
    } else {
        vec.push_back(0);
    }
}


vector<double> tally_counts(const Parameters* par, Community* community, int discard_days) {
    assert((int) par->runLength >= discard_days + OVERRUN);                 // number of sim days for data aggregation must be greater than OVERRUN + days discarded (from the beginning of sim)
    const size_t num_weeks = (par->runLength - discard_days - OVERRUN)/7;   // number of full weeks of data to be aggregated

    //vector<size_t> infected    = community->getNumNewlyInfected();
    vector<size_t> symptomatic = community->getNumNewlySymptomatic();
    vector<size_t> severe      = community->getNumNewlySevere();
    vector<size_t> dead        = community->getNumNewlyDead();

    // pair of num of primary infections starting this day, and mean num secondary infections they cause
    vector<pair<size_t, double>> R = community->getMeanNumSecondaryInfections();
    vector<size_t> Rt_incidence_tally(num_weeks, 0);

    vector<double> metrics(num_weeks*4, 0.0);
    for (size_t t = discard_days; t < discard_days + (7*num_weeks); ++t) {
        const size_t w = (t-discard_days)/7; // which reporting week are we in?
        metrics[w]                 += symptomatic[t];
        metrics[num_weeks + w]     += severe[t];
        metrics[2 * num_weeks + w] += dead[t];
        metrics[3 * num_weeks + w] += R[t].first > 0 ? R[t].first*R[t].second : 0;
        Rt_incidence_tally[w]      += R[t].first;
    }

    for (size_t w = 0; w < num_weeks; ++w) {
        metrics[3 * num_weeks + w] /= Rt_incidence_tally[w] > 0 ? Rt_incidence_tally[w] : 1.0;
    }
//metrics.resize(300);
    return metrics;
}

vector<double> calc_Rt_moving_average(vector<pair<size_t, double>> Rt_pairs, size_t window) {
    assert(window % 2 == 1);
    vector<double> Rt_ma(Rt_pairs.size(), 0.0);

    size_t halfwindow = (window - 1)/2;
    for (size_t i = 0; i < Rt_ma.size(); ++i) {
        size_t window_ct = 0;
        double infections = 0;
        for (size_t wi = halfwindow >= i ? 0 : i - halfwindow;
             wi < (i + halfwindow + 1 >= Rt_pairs.size() ? Rt_pairs.size() : i + halfwindow + 1);
             ++wi) {
            const size_t ct = Rt_pairs[wi].first;
            const double Rt = Rt_pairs[wi].second;
            window_ct += ct;
            infections += ct*Rt;
        }
        Rt_ma[i] = infections / window_ct;
//        cerr << "day, inf, ct, Rt_ma, Rt_count, Rt_val: " << i << " " << infections << " " << window_ct << " " << Rt_ma[i] << " " << Rt_pairs[i].first << " " << Rt_pairs[i].second << endl;
    }
    return Rt_ma;
}


#endif
