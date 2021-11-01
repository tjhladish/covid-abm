#include "Vac_Campaign.h"
#include "Community.h"
#include "Person.h"
#include "Location.h"

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

void Vac_Campaign::geographic_scheduling(int day, set<Person*> targetedPeople, Community* community) {
    if (day >= start_of_campaign[GEO_VACCINATION]) {
        for (Person* p : targetedPeople) {
            Location* center_loc = p->getHomeLoc();
            double radius = 0.001;                                  // distance in decimal degrees to add to center_loc to draw the capture area
            double northern_bound = center_loc->getY() + radius;
            double southern_bound = center_loc->getY() - radius;
            double western_bound  = center_loc->getX() + radius;
            double eastern_bound  = center_loc->getX() - radius;

            set<Location*> targetedLocs;
            for (Location* loc : community->getLocations()) {
                if (loc->getY() >= southern_bound and loc->getY() <= northern_bound and
                    loc->getX() >= eastern_bound  and loc->getX() <= western_bound) { targetedLocs.insert(loc); }
            }

            double coverage = 1.0;
            for (Location* loc : targetedLocs) {
                for (Person* p : loc->getPeople()) {
                    if (gsl_rng_uniform(RNG) < coverage and is_age_eligible_on(p->getAge(), day)
                        and not (p->hasBeenInfected() or p->isVaccinated())) { prioritize_vaccination(p); }
                }
            }
        }
    }
}

void Vac_Campaign::location_scheduling(int day, set<Person*> targetedPeople) {
    if (day >= start_of_campaign[LOCATION_VACCINATION]) {
        double coverage = 1.0;
        for (Person* tp : targetedPeople) {
            for (Person* p : tp->getHomeLoc()->getPeople()) {
                if (gsl_rng_uniform(RNG) < coverage and is_age_eligible_on(p->getAge(), day)
                    and not (p->hasBeenInfected() or p->isVaccinated())) { prioritize_vaccination(p); }
            }
            for (Person* p : tp->getDayLoc()->getPeople()) {
                if (gsl_rng_uniform(RNG) < coverage and is_age_eligible_on(p->getAge(), day)
                    and not (p->hasBeenInfected() or p->isVaccinated())) { prioritize_vaccination(p); }
            }
        }
    }
}
