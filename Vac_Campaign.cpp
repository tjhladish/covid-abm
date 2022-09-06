#include "Vac_Campaign.h"
#include "Community.h"
#include "Person.h"
#include "Location.h"
#include "Utility.h"
#include <map>

using covid::util::choose_k;

vector<Eligibility_Group*> Vac_Campaign::init_new_eligible_groups(int day) {
    // create empty eligibility group to add new urgent people to the queue
    vector<Eligibility_Group*> eligibles(_par->numVaccineDoses);
    for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
        Eligibility_Group* eg = new Eligibility_Group();
        eg->eligibility_day = (dose == 0) ? day : day + _par->vaccineDoseInterval.at(dose - 1);
        eligibles[dose] = eg;
    }

    return eligibles;
}

// TODO: will people be urgently queued for doses other than dose 1?
void Vac_Campaign::ring_scheduling(int day, vector<set<Person*, PerPtrComp>> tracedContacts) {
    if (day >= start_of_campaign[RING_VACCINATION]) {
        //schedule urgent vax
        // TODO add nuaince to...
        //  - the decision of whether someone should be sch based on time since last infection
        //  - how we decide/set which groups will be scheduled

        // create empty eligibility group to add new urgent people to the queue
        vector<Eligibility_Group*> urgents = init_new_eligible_groups(day);

        for (auto s : tracedContacts) {
            for (Person* p : s) {
                // Parameters must ensure than urgent_vax_dose_threshold <= numVaccineDoses
                // if (is_age_eligible_on(p->getAge(), day) and not (p->hasBeenInfected() or (p->getNumVaccinations() < _par->urgent_vax_dose_threshold))) {
                if (is_age_eligible_on(p->getAge(), day) and (p->getNumVaccinations() < _par->urgent_vax_dose_threshold)) {
                    // getNumVaccinations() used as an index will ensure this person is being scheduled for their next dose
                    // eg: if getNumVaccinations() returns 1, 1 as an index pushes this person for dose 2
                    urgents[p->getNumVaccinations()]->eligible_people[age_bin_lookup[p->getAge()]].push_back(p);
                }
            }
        }
        schedule_urgent_doses(urgents);
    }
    return;
}

void Vac_Campaign::geographic_scheduling(int day, vector<set<Person*, PerPtrComp>> targetedPeople, Community* community) {
    if (day >= start_of_campaign[GEO_VACCINATION]) {
        set<Person*, PerPtrComp> tracedCases = targetedPeople[0];
        set<Location*> tracedLocs;
        for (Person* p : tracedCases) { tracedLocs.insert(p->getHomeLoc()); }

        set<Location*> targetedLocs;
        const double radius = 0.001;    // distance in decimal degrees to add to center_loc to draw the capture area
        for (Location* center_loc : tracedLocs) {
            vector<Location*> locs_to_search;
            for (double y = -0.01; y <= 0.01; y += 0.01) {
                for (double x = -0.01; x <= 0.01; x += 0.01) {
                    vector<Location*> tmp = community->locsAtPixel(center_loc->getXPixel() + x, center_loc->getYPixel() + y);
                    locs_to_search.insert(locs_to_search.end(), tmp.begin(), tmp.end());
                }
            }
            //vector<Location*> locs_to_search = community->locsAtPixel(center_loc->getXPixel(), center_loc->getYPixel());

            const double northern_bound = center_loc->getY() + radius;
            const double southern_bound = center_loc->getY() - radius;
            const double eastern_bound  = center_loc->getX() + radius;
            const double western_bound  = center_loc->getX() - radius;

            for (Location* loc : locs_to_search) {
                if (loc->getY() > southern_bound and loc->getY() < northern_bound and
                    loc->getX() < eastern_bound  and loc->getX() > western_bound) { targetedLocs.insert(loc); }
            }
        }

        // create empty eligibility group to add new urgent people to the queue
        vector<Eligibility_Group*> urgents = init_new_eligible_groups(day);

        double coverage = 1.0;
        for (Location* loc : targetedLocs) {
            for (Person* p : loc->getPeople()) {
                if (gsl_rng_uniform(RNG) < coverage and is_age_eligible_on(p->getAge(), day)
                    and not (p->hasBeenInfected() or p->isVaccinated())) {
                    // getNumVaccinations() used as an index will ensure this person is being scheduled for their next dose
                    // eg: if getNumVaccinations() returns 1, 1 as an index pushes this person for dose 2
                    urgents[p->getNumVaccinations()]->eligible_people[age_bin_lookup[p->getAge()]].push_back(p);
                }
            }
        }
        schedule_urgent_doses(urgents);
    }
}

void Vac_Campaign::location_scheduling(int day, vector<set<Person*, PerPtrComp>> targetedPeople) {
    if (day >= start_of_campaign[LOCATION_VACCINATION]) {
        const double work_and_neighbr_coverage = 0.5;
        const size_t expected_classroom_size   = 30;

        set<Person*, PerPtrComp> tracedCases = targetedPeople[0];
        set<Person*, PerPtrComp> all_people_to_vax;

        // create empty eligibility group to add new urgent people to the queue
        vector<Eligibility_Group*> urgents = init_new_eligible_groups(day);

        for (Person* tp : tracedCases) {
            Location* homeloc = tp->getHomeLoc();
            for (Person* p : homeloc->getPeople()) { all_people_to_vax.insert(p); }

            vector<Person*> all_neighbors;
            for (Location* n : homeloc->getNeighbors()) {
                for (Person* p : n->getPeople()) { all_neighbors.push_back(p); }
            }
            const size_t num_neighbors_to_vax = all_neighbors.size() * work_and_neighbr_coverage;
            vector<Person*> neighbors_to_vax = choose_k(RNG, all_neighbors, num_neighbors_to_vax);
            for (Person* p : neighbors_to_vax) { all_people_to_vax.insert(p); }

            if (tp->getDayLoc()) {
                Location* dayloc = tp->getDayLoc();
                vector<Person*> dayloc_people_to_vax;

                if (dayloc->getType() == WORK) {
                    const size_t num_workers_to_vax = dayloc->getNumPeople() * work_and_neighbr_coverage;
                    vector<Person*> workplace_people = dayloc->getPeople();
                    dayloc_people_to_vax = choose_k(RNG, workplace_people, num_workers_to_vax);
                } else if (dayloc->getType() == SCHOOL) {
                    size_t num_classmates_to_vax = gsl_ran_poisson(RNG, expected_classroom_size);
                    num_classmates_to_vax = (num_classmates_to_vax > (size_t) dayloc->getNumPeople()) ? dayloc->getNumPeople() : num_classmates_to_vax;
                    vector<Person*> school_people = dayloc->getPeople();
                    dayloc_people_to_vax = choose_k(RNG, school_people, num_classmates_to_vax);
                }

                for (Person* p : dayloc_people_to_vax) { all_people_to_vax.insert(p); }
            }

            for (Person* p : all_people_to_vax) {
                if (is_age_eligible_on(p->getAge(), day) and not (p->hasBeenInfected() or p->isVaccinated())) {
                    // getNumVaccinations() used as an index will ensure this person is being scheduled for their next dose
                    // eg: if getNumVaccinations() returns 1, 1 as an index pushes this person for dose 2
                    urgents[p->getNumVaccinations()]->eligible_people[age_bin_lookup[p->getAge()]].push_back(p);
                }
            }
        }
        schedule_urgent_doses(urgents);
    }
}

void Vac_Campaign::grouped_risk_scheduling(int day, Community* community) {
    // this function only needs to be executed once at the start of the strategy
    // grouping and scheduling is performed here and then the groups will be added by specialty methods
    if (day == start_of_campaign[GROUPED_RISK_VACCINATION]) {
        map<int, double> grouped_risk;   // map of group to pair of risk sum and group size (can then calculate per capita risk)
        map<int, vector<Person*>> grouped_ppl;      // map of group to vector of people in that group
        ifstream iss(_par->riskGroupsFilename);

        if (!iss) { cerr << "ERROR: " << _par->riskGroupsFilename << " not found." << endl; }

        string buffer;
        istringstream line;
        int pid, group;

        // the input file is expected to have two columns: person ID and the group they belong to
        while ( getline(iss, buffer) ) {
            line.clear();
            line.str(buffer);

            // for each line in the file, calculate the individual's risk for their group and add
            if (line >> pid >> group) {
                Person* p = community->getPersonByID(pid);
                grouped_ppl[group].push_back(p);

                // risk of severe outcomes given infection
                double risk = _par->probSeriousOutcome.at(SEVERE)[p->hasComorbidity()][p->getAge()] * _par->pathogenicityByAge[p->getAge()];
                if (not grouped_risk.count(group)) {
                    grouped_risk[group]  = 0.0;
                }

                grouped_risk[group] += risk;
            }
        }
        iss.close();

        // special comparator to sort the groups by per capita risk
        struct group_risk_comp {
            bool operator() (pair<int, double> const& lhs, pair<int, double> const& rhs) const {
                return lhs.second > rhs.second;
            }
        };

        // for each group of people, calculate the per capita risk of the group and sort the groups by per capita risk
        set< pair<int, double>, group_risk_comp > ordered_groups_by_per_capita_risk;
        for (const auto [group, totRisk] : grouped_risk) {
            pair<int, double> tmp(group, (double) totRisk/grouped_ppl[group].size());
            ordered_groups_by_per_capita_risk.insert(tmp);
        }

        // for each group create a new eligibility group for the group that can be dumped into the vac_campaign urgent pool
        _grouped_risk_deque.clear();
        for (pair<int, double> group : ordered_groups_by_per_capita_risk) {
            vector<Eligibility_Group*> urgents = init_new_eligible_groups(day);
            for (Person* p : grouped_ppl[group.first]) {
                // Parameters must ensure than urgent_vax_dose_threshold <= numVaccineDoses
                // if (is_age_eligible_on(p->getAge(), day) and not (p->hasBeenInfected() or (p->getNumVaccinations() < _par->urgent_vax_dose_threshold))) {
                if (is_age_eligible_on(p->getAge(), day) and (p->getNumVaccinations() < _par->urgent_vax_dose_threshold)) {
                    _sch_risk_groups[group.first].push_back(p);
                    // getNumVaccinations() used as an index will ensure this person is being scheduled for their next dose
                    // eg: if getNumVaccinations() returns 1, 1 as an index pushes this person for dose 2
                    urgents[p->getNumVaccinations()]->eligible_people[age_bin_lookup[p->getAge()]].push_back(p);
                }
            };
            pair<int, vector<Eligibility_Group*>> tmp_risk_group(group.first, urgents);
            _grouped_risk_deque.push_back(tmp_risk_group);
        }

        // add first group to pool
        for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
            for (int bin : unique_age_bins) {
                _insert_eligible_people(_grouped_risk_deque.front().second[dose], potential_urg_vaccinees, dose, bin);
            }
            delete _grouped_risk_deque.front().second[dose];
        }
        _current_risk_group = _grouped_risk_deque.front().first;
        _grouped_risk_deque.pop_front();
    }
}

void Vac_Campaign::risk_scheduling(int day, Community* community, double hosp_risk_threshold) {
    // this function only needs to be executed once at the start of the strategy
    if (day == start_of_campaign[RISK_VACCINATION]) {
        vector<Eligibility_Group*> urgents = init_new_eligible_groups(day);

        for (Person* p : community->getPeople()) {
            // risk of severe outcomes given infection
            double risk = _par->probSeriousOutcome.at(SEVERE)[p->hasComorbidity()][p->getAge()] * _par->pathogenicityByAge[p->getAge()];
            if (risk >= hosp_risk_threshold and is_age_eligible_on(p->getAge(), day) and (p->getNumVaccinations() < _par->urgent_vax_dose_threshold)) {
                // getNumVaccinations() used as an index will ensure this person is being scheduled for their next dose
                // eg: if getNumVaccinations() returns 1, 1 as an index pushes this person for dose 2
                urgents[p->getNumVaccinations()]->eligible_people[age_bin_lookup[p->getAge()]].push_back(p);
            }
        }

        schedule_urgent_doses(urgents);
    }
}

void Vac_Campaign::reactive_strategy(int day, vector<set<Person*, PerPtrComp>> targetedPeople, Community* community) {
    // general handler for which campaign should be called from Community::tick()
    switch (reactive_vac_strategy) {
        case RING_VACCINATION:         { ring_scheduling(day, targetedPeople); break; }
        case GEO_VACCINATION:          { geographic_scheduling(day, targetedPeople, community); break; }
        case LOCATION_VACCINATION:     { location_scheduling(day, targetedPeople); break; }
        case GROUPED_RISK_VACCINATION: { grouped_risk_scheduling(day, community); break; }
        case RISK_VACCINATION:         { risk_scheduling(day, community); break; }

        case GENERAL_CAMPAIGN: break;
        case NUM_OF_VAC_CAMPAIGN_TYPES: break;
        default: break;
    }
}

void Vac_Campaign::init_eligibility_queue(const Community* community) {
    std_eligibility_queue.clear();
    std_eligibility_queue.resize(_par->numVaccineDoses);

    urg_eligibility_queue.clear();
    urg_eligibility_queue.resize(_par->numVaccineDoses);

    Eligibility_Group* first_eg = new Eligibility_Group();
    first_eg->eligibility_day = 0;

    // everyone in pop becomes eligible for first dose on day 0
    for (Person* p : community->getPeople()) {
        first_eg->eligible_people[age_bin_lookup[p->getAge()]].push_back(p);
    }
    std_eligibility_queue[0].push(first_eg);
}

// generates and saves comprehensive, mutually exclusive age bins for the entire population
void Vac_Campaign::generate_age_bins(Community* community, std::set<int> unique_bin_mins, std::set<int> unique_bin_maxs) {
    std::vector<int> mins(unique_bin_mins.begin(), unique_bin_mins.end());
    std::vector<int> maxs(unique_bin_maxs.begin(), unique_bin_maxs.end());

    std::set<int> ages_not_binned;
    assert(mins.size() == maxs.size());

    // for every possible age, check if that age is covered by an age bin from the input file
    // after, only ages needing to be binned will remain
    for (int age = 0; age < NUM_AGE_CLASSES; ++age) {
        ages_not_binned.insert(age);
        for (int i = 0; i < (int)mins.size(); ++i) {
            if ((age >= mins[i]) and (age <= maxs[i])) {
                ages_not_binned.erase(age);
                break;
            }
        }
    }

    // for remaining ages needing to be binnned, find their bins
    // will check each value for contiguity
    std::vector<int> ages_to_be_binned(ages_not_binned.begin(), ages_not_binned.end());
    for (int i = 0; i < (int)ages_to_be_binned.size(); ++i) {
        if (i == 0) { mins.push_back(ages_to_be_binned[i]); continue; }

        if (i == ((int)ages_to_be_binned.size() - 1)) { maxs.push_back(ages_to_be_binned[i]); }

        if (ages_to_be_binned[i] == (ages_to_be_binned[i - 1] + 1)) {
            continue;
        } else {
            mins.push_back(ages_to_be_binned[i]);
            maxs.push_back(ages_to_be_binned[i - 1]);
        }
    }

    assert(mins.size() == maxs.size());
    for (int i = 0; i < (int)mins.size(); ++i) { unique_age_bin_pops[mins[i]] = 0; }

    // for each age, save the bin min associated with it for fast lookup
    // for each bin, aggregate and save the bin population
    age_bin_lookup = std::vector<int>(NUM_AGE_CLASSES);
    for (int age = 0; age < NUM_AGE_CLASSES; ++age) {
        for (int i = 0; i < (int)mins.size(); ++i) {
            if ((age >= mins[i]) and (age <= maxs[i])) {
                age_bin_lookup[age] = mins[i];
                unique_age_bin_pops[mins[i]] += community->getAgeCohort(age).size();
                break;
            }
        }
    }

    unique_age_bins = mins;
}

// called when doses should be aggregated regardless of dose or age bin
// all dose, bin combinations for a given day will point to the same int in _doses and share those doses
Dose_Ptrs Vac_Campaign::_dose_pool(Dose_Vals doses_in) {
    Dose_Ptrs doses_available(_par->runLength, std::vector< std::map<int, int*> >(_par->numVaccineDoses));
    for (int day = 0; day < (int) _par->runLength; ++day) {
        _doses.push_back(0);
        int* daily_pooled_dose_ptr = &(_doses.back());
        for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
            for (int bin : unique_age_bins) {
                *daily_pooled_dose_ptr += doses_in.at(day).at(dose).at(bin);
                doses_available.at(day).at(dose)[bin] = daily_pooled_dose_ptr;
            }
        }
    }
    return doses_available;
}

// doses pooled across std and urg
void Vac_Campaign::_dose_pool_all(vector<Dose_Ptrs> &doses_out, /*Dose_Ptrs &std_out, */Dose_Vals urg_in, Dose_Vals std_in) {
    Dose_Ptrs doses_available(_par->runLength, std::vector< std::map<int, int*> >(_par->numVaccineDoses));
    for (int day = 0; day < (int) _par->runLength; ++day) {
        _doses.push_back(0);
        int* daily_pooled_dose_ptr = &(_doses.back());
        for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
            for (int bin : unique_age_bins) {
                *daily_pooled_dose_ptr += urg_in.at(day).at(dose).at(bin) + std_in.at(day).at(dose).at(bin);
                doses_available.at(day).at(dose)[bin] = daily_pooled_dose_ptr;
            }
        }
    }
    doses_out[URGENT_ALLOCATION]   = doses_available;
    doses_out[STANDARD_ALLOCATION] = doses_available;
}

// stores doses as read in the input file
// each day, dose, bin combinations will point to a unique int in _doses representing that combination's doses available
Dose_Ptrs Vac_Campaign::_dose_store(Dose_Vals doses_in) {
    Dose_Ptrs doses_available(_par->runLength, std::vector< std::map<int, int*> >(_par->numVaccineDoses));
    for (int day = 0; day < (int) _par->runLength; ++day) {
        for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
            for (int bin : unique_age_bins) {
                _doses.push_back(doses_in.at(day).at(dose).at(bin));
                doses_available.at(day).at(dose)[bin] = &(_doses.back());
            }
        }
    }
    return doses_available;
}

void Vac_Campaign::_clear_and_resize_doses() {
    // max number of _doses elements is the total possible number of day, dose, bin combinations for standard and urgent (hence the x2)
    const size_t max_num_elements = _par->runLength * _par->numVaccineDoses * unique_age_bins.size() * 2;
    // the vector must be of size 0 but have the capacity for max_num_elements to prevent reallocation of the vector (nullifying all former pointers to it)
    _doses = std::vector<int>(0);
    _doses.reserve(max_num_elements);
}

void Vac_Campaign::init_doses_available(Dose_Vals urg_in, Dose_Vals std_in) {
    _clear_and_resize_doses();
    if (pool_urg_doses and not pool_std_doses) {
        doses_available[URGENT_ALLOCATION]   = _dose_pool(urg_in);   // pool urgent doses regardless of age or dose
        doses_available[STANDARD_ALLOCATION] = _dose_store(std_in);  // store standard doses normally
    } else if (pool_std_doses and not pool_urg_doses) {
        doses_available[URGENT_ALLOCATION]   = _dose_store(urg_in);  // store urgent doses normally
        doses_available[STANDARD_ALLOCATION] = _dose_pool(std_in);   // pool standard doses regardless of age or dose
    } else if (pool_urg_doses and pool_std_doses) {
        // pool all doses regardless of age or dose
        doses_available[URGENT_ALLOCATION]   = _dose_pool(urg_in);
        doses_available[STANDARD_ALLOCATION] = _dose_pool(std_in);
    } else if (pool_all_doses) {
        // need a new function to pool across std and urg
        _dose_pool_all(doses_available, urg_in, std_in);
    } else {
        // store all doses normally
        doses_available[URGENT_ALLOCATION]   = _dose_store(urg_in);
        doses_available[STANDARD_ALLOCATION] = _dose_store(std_in);
    }
}

// called from the Community copy ctor to handle the construction of a new _doses to ensure it does not move and invalidate any ptrs
void Vac_Campaign::copy_doses_available(Vac_Campaign* ovc) {
    _clear_and_resize_doses();
    std::map<int*, int*> dose_ptr_map;
    for (int day = 0; day < (int) _par->runLength; ++day) {
        for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
            for (int bin : unique_age_bins) {
                int* other_urg_ptr = ovc->doses_available[URGENT_ALLOCATION][day][dose][bin];
                int* other_std_ptr = ovc->doses_available[STANDARD_ALLOCATION][day][dose][bin];

                if (dose_ptr_map.count(other_urg_ptr)) { // TODO - this can be refactored to be simpler
                    // ptr has already been seen; use the same new ptr to maintain pooled relationships
                    doses_available[URGENT_ALLOCATION][day][dose][bin] = dose_ptr_map[other_urg_ptr];
                } else {
                    // this is a new ptr; push back value and create new ptr
                    _doses.push_back(*other_urg_ptr);
                    doses_available[URGENT_ALLOCATION][day][dose][bin] = &(_doses.back());
                    dose_ptr_map[other_urg_ptr] = &(_doses.back());
                }

                if (dose_ptr_map.count(other_std_ptr)) {
                    // ptr has already been seen; use the same new ptr to maintain pooled relationships
                    doses_available[STANDARD_ALLOCATION][day][dose][bin] = dose_ptr_map[other_std_ptr];
                } else {
                    // this is a new ptr; push back value and create new ptr
                    _doses.push_back(*other_std_ptr);
                    doses_available[STANDARD_ALLOCATION][day][dose][bin] = &(_doses.back());
                    dose_ptr_map[other_std_ptr] = &(_doses.back());
                }
            }
        }
    }
}

Vac_Campaign* Vac_Campaign::quick_cache() {
    Vac_Campaign* vc       = new Vac_Campaign(*this);
    Eligibility_Q other_sq = this->get_std_eligibility_queue();
    Eligibility_Q other_uq = this->get_urg_eligibility_queue();

    Eligibility_Q sq = vc->get_std_eligibility_queue();
    Eligibility_Q uq = vc->get_urg_eligibility_queue();
    sq.clear(); sq.resize(_par->numVaccineDoses);
    uq.clear(); uq.resize(_par->numVaccineDoses);

    for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
        while (other_sq[dose].size()) {
            Eligibility_Group* other_eg = other_sq[dose].top();
            Eligibility_Group* eg = new Eligibility_Group();
            eg->eligibility_day = other_eg->eligibility_day;
            for (int bin : vc->get_unique_age_bins()) {
                for (Person* p : other_eg->eligible_people[bin]) {
                    eg->eligible_people[bin].push_back(p);
                }
            }
            sq[dose].push(eg);
            other_sq[dose].pop();
        }

        while (other_uq[dose].size()) {
            Eligibility_Group* other_eg = other_uq[dose].top();
            Eligibility_Group* eg = new Eligibility_Group();
            eg->eligibility_day = other_eg->eligibility_day;
            for (int bin : vc->get_unique_age_bins()) {
                for (Person* p : other_eg->eligible_people[bin]) {
                    eg->eligible_people[bin].push_back(p);
                }
            }
            uq[dose].push(eg);
            other_uq[dose].pop();
        }
    }

    vc->copy_doses_available(this);
    vc->set_std_eligibility_queue(sq);
    vc->set_urg_eligibility_queue(uq);

    return vc;
}
