#include "Vac_Campaign.h"
#include "Community.h"
#include "Person.h"
#include "Location.h"
#include "Utility.h"

using covid::util::choose_k;

// void Vac_Campaign::ring_scheduling(int day, vector< set<Person*> > tracedContacts) {
//     if (day >= start_of_campaign[RING_VACCINATION]) {
//         //schedule urgent vax
//         // TODO add nuaince to...
//         //  - the decision of whether someone should be sch based on time since last infection
//         //  - how we decide/set which groups will be scheduled
//         for (auto s : tracedContacts) {
//             for (Person* p : s) {
//                 if (is_age_eligible_on(p->getAge(), day) and not (p->hasBeenInfected() or p->isVaccinated())) { prioritize_vaccination(p); }
//             }
//         }
//     }
//     return;
// }
//
// void Vac_Campaign::geographic_scheduling(int day, vector< set<Person*> > targetedPeople, Community* community) {
//     if (day >= start_of_campaign[GEO_VACCINATION]) {
//         set<Person*> tracedCases = targetedPeople[0];
//         set<Location*> tracedLocs;
//         for (Person* p : tracedCases) { tracedLocs.insert(p->getHomeLoc()); }
//
//         set<Location*> targetedLocs;
//         const double radius = 0.001;                                  // distance in decimal degrees to add to center_loc to draw the capture area
//         for (Location* center_loc : tracedLocs) {
//             vector<Location*> locs_to_search;
//             for (double y = -0.01; y <= 0.01; y += 0.01) {
//                 for (double x = -0.01; x <= 0.01; x += 0.01) {
//                     vector<Location*> tmp = community->locsAtPixel(center_loc->getXPixel() + x, center_loc->getYPixel() + y);
//                     locs_to_search.insert(locs_to_search.end(), tmp.begin(), tmp.end());
//                 }
//             }
//             //vector<Location*> locs_to_search = community->locsAtPixel(center_loc->getXPixel(), center_loc->getYPixel());
//
//             const double northern_bound = center_loc->getY() + radius;
//             const double southern_bound = center_loc->getY() - radius;
//             const double eastern_bound  = center_loc->getX() + radius;
//             const double western_bound  = center_loc->getX() - radius;
//
//             for (Location* loc : locs_to_search) {
//                 if (loc->getY() > southern_bound and loc->getY() < northern_bound and
//                     loc->getX() < eastern_bound  and loc->getX() > western_bound) { targetedLocs.insert(loc); }
//             }
//         }
//         //for (Person* p : targetedPeople) {
//         //for (Person* p : tracedCases) {
//         //    Location* center_loc = p->getHomeLoc();
//         //    double radius = 0.001;                                  // distance in decimal degrees to add to center_loc to draw the capture area
//         //    double northern_bound = center_loc->getY() + radius;
//         //    double southern_bound = center_loc->getY() - radius;
//         //    double western_bound  = center_loc->getX() + radius;
//         //    double eastern_bound  = center_loc->getX() - radius;
//
//         //    for (Location* loc : community->getLocations()) {
//         //        if (loc->getY() >= southern_bound and loc->getY() <= northern_bound and
//         //            loc->getX() >= eastern_bound  and loc->getX() <= western_bound) { targetedLocs.insert(loc); }
//         //    }
//         //}
//
//         double coverage = 1.0;
//         for (Location* loc : targetedLocs) {
//             for (Person* p : loc->getPeople()) {
//                 if (gsl_rng_uniform(RNG) < coverage and is_age_eligible_on(p->getAge(), day)
//                     and not (p->hasBeenInfected() or p->isVaccinated())) { prioritize_vaccination(p); }
//             }
//         }
//     }
// }
//
// void Vac_Campaign::location_scheduling(int day, vector< set<Person*> > targetedPeople) {
//     if (day >= start_of_campaign[LOCATION_VACCINATION]) {
//         const double work_and_neighbr_coverage = 0.5;
//         const size_t expected_classroom_size   = 30;
//
//         set<Person*> tracedCases = targetedPeople[0];
//         set<Person*> all_people_to_vax;
//         //for (Person* tp : targetedPeople) {
//         for (Person* tp : tracedCases) {
//             Location* homeloc = tp->getHomeLoc();
//             for (Person* p : homeloc->getPeople()) { all_people_to_vax.insert(p); }
//
//             vector<Person*> all_neighbors;
//             for (Location* n : homeloc->getNeighbors()) {
//                 for (Person* p : n->getPeople()) { all_neighbors.push_back(p); }
//             }
//             const size_t num_neighbors_to_vax = all_neighbors.size() * work_and_neighbr_coverage;
//             vector<Person*> neighbors_to_vax = choose_k(RNG, all_neighbors, num_neighbors_to_vax);
//             for (Person* p : neighbors_to_vax) { all_people_to_vax.insert(p); }
//
//             if (tp->getDayLoc()) {
//                 Location* dayloc = tp->getDayLoc();
//                 vector<Person*> dayloc_people_to_vax;
//
//                 if (dayloc->getType() == WORK) {
//                     const size_t num_workers_to_vax = dayloc->getNumPeople() * work_and_neighbr_coverage;
//                     vector<Person*> workplace_people = dayloc->getPeople();
//                     dayloc_people_to_vax = choose_k(RNG, workplace_people, num_workers_to_vax);
//                 } else if (dayloc->getType() == SCHOOL) {
//                     size_t num_classmates_to_vax = gsl_ran_poisson(RNG, expected_classroom_size);
//                     num_classmates_to_vax = (num_classmates_to_vax > (size_t) dayloc->getNumPeople()) ? dayloc->getNumPeople() : num_classmates_to_vax;
//                     vector<Person*> school_people = dayloc->getPeople();
//                     dayloc_people_to_vax = choose_k(RNG, school_people, num_classmates_to_vax);
//                 }
//
//                 for (Person* p : dayloc_people_to_vax) { all_people_to_vax.insert(p); }
//             }
//
//             for (Person* p : all_people_to_vax) {
//                 if (is_age_eligible_on(p->getAge(), day) and not (p->hasBeenInfected() or p->isVaccinated())) { prioritize_vaccination(p); }
//             }
//         }
//     }
// }
//
// void Vac_Campaign::reactive_strategy(int day, vector< set<Person*> > targetedPeople, Community* community) {
//     // general handler for which campgain shoudl be called from Community::tick()
//     switch (reactive_vac_strategy) {
//         case RING_VACCINATION:     { ring_scheduling(day, targetedPeople); break; }
//         case GEO_VACCINATION:      { geographic_scheduling(day, targetedPeople, community); break; }
//         case LOCATION_VACCINATION: { location_scheduling(day, targetedPeople); break; }
//
//         case GENERAL_CAMPAIGN: break;
//         case NUM_OF_VAC_CAMPAIGN_TYPES: break;
//         default: break;
//     }
// }

void Vac_Campaign::init_eligibility_queue(const Community* community) {
    eligibility_queue.clear();
    eligibility_queue.resize(_par->numVaccineDoses);
    Eligibility_Group* first_eg = new Eligibility_Group();
    first_eg->eligibility_day = 0;

    // everyone in pop becomes eligible for first dose on day 0
    for (Person* p : community->getPeople()) {
        first_eg->eligible_people[age_bin_lookup[p->getAge()]].push_back(p);
    }
    eligibility_queue[0].push(first_eg);
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
