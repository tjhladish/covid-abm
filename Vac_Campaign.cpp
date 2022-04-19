#include "Vac_Campaign.h"
#include "Community.h"
#include "Person.h"
#include "Location.h"
#include "Utility.h"

using covid::util::choose_k;

vector<Eligibility_Group*> Vac_Campaign::init_new_eligible_groups(int day) {
    // create empty eligibility group to add new urgent people to the queue
    vector<Eligibility_Group*> eligibles(_par->numVaccineDoses);
    for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
        Eligibility_Group* eg = new Eligibility_Group();
        eg->eligibility_day = (dose == 0) ? day + 1 : day + _par->vaccineDoseInterval.at(dose - 1);
        eligibles[dose] = eg;
    }

    return eligibles;
}

// TODO: will people be urgently queued for doses other than dose 1?
void Vac_Campaign::ring_scheduling(int day, vector< set<Person*> > tracedContacts) {
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

void Vac_Campaign::geographic_scheduling(int day, vector< set<Person*> > targetedPeople, Community* community) {
    if (day >= start_of_campaign[GEO_VACCINATION]) {
        set<Person*> tracedCases = targetedPeople[0];
        set<Location*> tracedLocs;
        for (Person* p : tracedCases) { tracedLocs.insert(p->getHomeLoc()); }

        set<Location*> targetedLocs;
        const double radius = 0.001;                                  // distance in decimal degrees to add to center_loc to draw the capture area
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

void Vac_Campaign::location_scheduling(int day, vector< set<Person*> > targetedPeople) {
    if (day >= start_of_campaign[LOCATION_VACCINATION]) {
        const double work_and_neighbr_coverage = 0.5;
        const size_t expected_classroom_size   = 30;

        set<Person*> tracedCases = targetedPeople[0];
        set<Person*> all_people_to_vax;

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

void Vac_Campaign::reactive_strategy(int day, vector< set<Person*> > targetedPeople, Community* community) {
    // general handler for which campgain shoudl be called from Community::tick()
    switch (reactive_vac_strategy) {
        case RING_VACCINATION:     { ring_scheduling(day, targetedPeople); break; }
        case GEO_VACCINATION:      { geographic_scheduling(day, targetedPeople, community); break; }
        case LOCATION_VACCINATION: { location_scheduling(day, targetedPeople); break; }

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
    std::vector ages_to_be_binned(ages_not_binned.begin(), ages_not_binned.end());
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

dosePtrs Vac_Campaign::_dose_pool(doseVals doses_in) {
    dosePtrs doses_available(_par->runLength, std::vector< std::map<int, int*> >(_par->numVaccineDoses));
    for (int day = 0; day < (int) _par->runLength; ++day) {
        _doses.push_back(0);
        int* daily_pooled_dose_ptr = &_doses.back();
        for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
            for (int bin : unique_age_bins) {
                *daily_pooled_dose_ptr += doses_in.at(day).at(dose).at(bin);
                doses_available.at(day).at(dose)[bin] = daily_pooled_dose_ptr;
            }
        }
    }
    return doses_available;
}

dosePtrs Vac_Campaign::_dose_store(doseVals doses_in) {
    dosePtrs doses_available(_par->runLength, std::vector< std::map<int, int*> >(_par->numVaccineDoses));
    for (int day = 0; day < (int) _par->runLength; ++day) {
        for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
            for (int bin : unique_age_bins) {
                _doses.push_back(0);
                int* dose_ptr = &_doses.back();
                *dose_ptr = doses_in.at(day).at(dose).at(bin);
                doses_available.at(day).at(dose)[bin] = dose_ptr;
            }
        }
    }
    return doses_available;
}

void Vac_Campaign::init_doses_available(doseVals urg_in, doseVals std_in) {
    // max number of _doses elements is the total possible number of day, dose, bin combinations for standard and urgent (hence the x2)
    const size_t max_num_elements = _par->runLength * _par->numVaccineDoses * unique_age_bins.size() * 2;
    _doses = std::vector<int>(0);
    _doses.reserve(max_num_elements);

    if (pool_urg_doses) {
        urg_doses_available = _dose_pool(urg_in);   // pool urgent doses regardless of age or dose
        std_doses_available = _dose_store(std_in);  // store standard doses normally
    } else if (pool_std_doses) {
        urg_doses_available = _dose_store(urg_in);  // store urgent doses normally
        std_doses_available = _dose_pool(std_in);   // pool standard doses regardless of age or dose
    } else if (pool_all_doses) {
        // pool all doses regardless of age or dose
        urg_doses_available = _dose_pool(urg_in);
        std_doses_available = _dose_pool(std_in);
    } else {
        // store all doses normally
        urg_doses_available = _dose_store(urg_in);
        std_doses_available = _dose_store(std_in);
    }
}
