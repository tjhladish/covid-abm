#include "Vac_Campaign.h"

void Vac_Campaign::ring_scheduling(int day, vector< set<Person*> > tracedContacts) {
    if(day >= start_of_campaign[RING_VACCINATION]){
        //schedule urgent vax
        // TODO add nuaince to...
        //  - the decision of whether someone should be sch based on time since last infection
        //  - how we decide/set which groups will be scheduled
        for(auto s : tracedContacts) {
            for(Person* p : s) {
                if(is_age_eligible_on(p->getAge(), day) and not (p->hasBeenInfected() or p->isVaccinated())) { prioritize_vaccination(p); }
            }
        }
    }
    return;
}
