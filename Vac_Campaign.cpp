#include "Vac_Campaign.h"
#include "Community.h"
#include "Person.h"
#include "Location.h"
#include "Utility.h"

using covid::util::choose_k;

void Vac_Campaign::ring_scheduling(int day, vector< set<Person*> > tracedContacts) {
    if (day >= start_of_campaign[RING_VACCINATION]) {
        //schedule urgent vax
        // TODO add nuaince to...
        //  - the decision of whether someone should be sch based on time since last infection
        //  - how we decide/set which groups will be scheduled
        for (auto s : tracedContacts) {
            for (Person* p : s) {
                if (is_age_eligible_on(p->getAge(), day) and not (p->hasBeenInfected() or p->isVaccinated())) { prioritize_vaccination(p); }
            }
        }
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
                const double northern_bound = center_loc->getY() + radius;
                const double southern_bound = center_loc->getY() - radius;
                const double eastern_bound  = center_loc->getX() + radius;
                const double western_bound  = center_loc->getX() - radius;

                vector<Location*> locs_to_search = community->locsAtPixel(center_loc->getXPixel(), center_loc->getYPixel());

            for (Location* loc : locs_to_search) {
                if (loc->getY() > southern_bound and loc->getY() < northern_bound and
                    loc->getX() < eastern_bound  and loc->getX() > western_bound) { targetedLocs.insert(loc); }
            }
        }
        //for (Person* p : targetedPeople) {
        //for (Person* p : tracedCases) {
        //    Location* center_loc = p->getHomeLoc();
        //    double radius = 0.001;                                  // distance in decimal degrees to add to center_loc to draw the capture area
        //    double northern_bound = center_loc->getY() + radius;
        //    double southern_bound = center_loc->getY() - radius;
        //    double western_bound  = center_loc->getX() + radius;
        //    double eastern_bound  = center_loc->getX() - radius;

        //    for (Location* loc : community->getLocations()) {
        //        if (loc->getY() >= southern_bound and loc->getY() <= northern_bound and
        //            loc->getX() >= eastern_bound  and loc->getX() <= western_bound) { targetedLocs.insert(loc); }
        //    }
        //}

        double coverage = 1.0;
        for (Location* loc : targetedLocs) {
            for (Person* p : loc->getPeople()) {
                if (gsl_rng_uniform(RNG) < coverage and is_age_eligible_on(p->getAge(), day)
                    and not (p->hasBeenInfected() or p->isVaccinated())) { prioritize_vaccination(p); }
            }
        }
    }
}

void Vac_Campaign::location_scheduling(int day, vector< set<Person*> > targetedPeople) {
    if (day >= start_of_campaign[LOCATION_VACCINATION]) {
        const double work_and_neighbr_coverage = 0.5;
        const size_t expected_classroom_size   = 30;

        set<Person*> tracedCases = targetedPeople[0];
        set<Person*> all_people_to_vax;
        //for (Person* tp : targetedPeople) {
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
                if (is_age_eligible_on(p->getAge(), day) and not (p->hasBeenInfected() or p->isVaccinated())) { prioritize_vaccination(p); }
            }
        }
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
