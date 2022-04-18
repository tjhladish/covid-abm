// Vac_Campaign.h
#ifndef __VAC_CAMPAIGN_H
#define __VAC_CAMPAIGN_H

#include "Person.h"
#include <deque>
#include <map>
#include <queue>
#include <set>

class Person;

enum VaccinationQueueType {
    URGENT_QUEUE,
    STANDARD_QUEUE,
    REVACCINATE_QUEUE,
    NUM_OF_VACCINATION_QUEUE_TYPES
};

enum VaccineAllocationType {
    URGENT_ALLOCATION,
    STANDARD_ALLOCATION,
    NUM_OF_VACCINE_ALLOCATION_TYPES
};


//enum ReactiveVaccinationStrategyType {
//    RING_VACCINATION,
//    GEOGRAPHIC_CLUSTER_VACCINATION,
//    NUM_OF_REACTIVE_VAC_STRATEGY_TYPES
//};

/*
 * forwardTracingGroups refers to the number of groups of people that will be contact traced
 *  - 1: reported cases for that day
 *  - 2: contacts (or primary contacts, level 1 contacts) of reported cases
 *  - 3: contacts of contacts (or secondary contacts, level 2 contacts)
 */
enum ContactTracingGroupType {
    REPORTED_CASES,
    PRIMARY_CONTACTS,
    SECONDARY_CONTACTS,
    NUM_OF_CONTACT_TRACING_GROUP_TYPES
};

enum VacCampaignType {
    GENERAL_CAMPAIGN,
    RING_VACCINATION,
    GEO_VACCINATION,
    LOCATION_VACCINATION,
    NUM_OF_VAC_CAMPAIGN_TYPES
};

class Vac_Campaign;
class Community;

class Eligibility_Group {
  public:
    Eligibility_Group() = default;
    Eligibility_Group(const Eligibility_Group& o) {
        eligibility_day = o.eligibility_day;
        eligible_people = o.eligible_people;
    }

    int eligibility_day;                                        // when these people become eligible for vaccination
    std::map<int, std::vector<Person*>> eligible_people;        // people to become eligible index by [age bin]

    // comparator to use to rank Eligibility_Groups by eligibility_day in a priority queue
    struct comparator { bool operator()(const Eligibility_Group* A, const Eligibility_Group* B) const { return A->eligibility_day > B->eligibility_day; } };
};

class Vaccinee {
    friend Vac_Campaign;
    public:
        Vaccinee() : person(nullptr), status(NUM_OF_VACCINATION_QUEUE_TYPES), dose_source(NUM_OF_VACCINE_ALLOCATION_TYPES) {}
        Vaccinee(Person* p, VaccinationQueueType stat, VaccineAllocationType source) : person(p), status(stat), dose_source(source) {}

        int getNumVaccinations () { return person->getNumVaccinations(); }
        Person* get_person() { return person; }
        VaccinationQueueType get_status() { return status; }
        VaccineAllocationType get_dose_source() { return dose_source; }

    private:
        Person* person;
        // We allow the possibility that the dose source won't match the vaccinee status because we may want to
        // consider idealized scenarios where doses can be rapidly reallocated between vaccine distribution
        // programs.
        VaccinationQueueType status;       // was this vaccinee standard (planned, proactive) or urgent (reactive, traced, etc.)
        VaccineAllocationType dose_source; // was this dose allocated to the urgent or standard queue

        bool vaccinate (int day) { return person->vaccinate(day); } // true if person was vaccinatable
};

class Vac_Campaign {
    public:
        Vac_Campaign(const Parameters* par) {
            _par = par;

            prioritize_first_doses = false;
            flexible_queue_allocation = false;
            unlim_urgent_doses = false;
            reactive_vac_strategy = NUM_OF_VAC_CAMPAIGN_TYPES;
            reactive_vac_dose_allocation = 0.0;

            start_of_campaign = vector<int>(NUM_OF_VAC_CAMPAIGN_TYPES, 0);
            end_of_campaign = vector<int>(NUM_OF_VAC_CAMPAIGN_TYPES, 0);

            orig_doses_available = std::vector< std::vector< std::map<int, int> > >(par->runLength, std::vector< std::map<int, int> >(par->numVaccineDoses));
            doses_available = std::vector< std::vector< std::map<int, int> > >(par->runLength, std::vector< std::map<int, int> >(par->numVaccineDoses));
            urg_doses_available = std::vector< std::vector< std::map<int, int> > >(par->runLength, std::vector< std::map<int, int> >(par->numVaccineDoses));

            doses_used = std::vector< std::vector< std::map<int, int> > >(par->runLength, std::vector< std::map<int, int> >(par->numVaccineDoses));
            urg_doses_used = std::vector< std::vector< std::map<int, int> > >(par->runLength, std::vector< std::map<int, int> >(par->numVaccineDoses));

            potential_vaccinees = std::vector< std::map<int, std::vector<Person*> > >(par->numVaccineDoses);
            urgent_vaccinees = std::vector< std::map<int, std::vector<Person*> > >(par->numVaccineDoses);

            eligibility_queue.clear();
            eligibility_queue.resize(_par->numVaccineDoses);

            urgent_queue.clear();
            urgent_queue.resize(_par->numVaccineDoses);
        }

        Vac_Campaign(const Vac_Campaign& o) {
            doses_available = o.doses_available;
            urg_doses_available = o.urg_doses_available;
            orig_doses_available = o.orig_doses_available;

            doses_used = o.doses_used;
            urg_doses_used = o.urg_doses_used;

            potential_vaccinees = o.potential_vaccinees;
            urgent_vaccinees = o.urgent_vaccinees;

            eligibility_queue = o.eligibility_queue;
            urgent_queue = o.urgent_queue;

            age_bin_lookup = o.age_bin_lookup;
            unique_age_bins = o.unique_age_bins;
            unique_age_bin_pops = o.unique_age_bin_pops;

            _par = o._par;
            _VAX_RNG = o._VAX_RNG;

            start_of_campaign = o.start_of_campaign;
            end_of_campaign = o.end_of_campaign;

            prioritize_first_doses = o.prioritize_first_doses;
            flexible_queue_allocation = o.flexible_queue_allocation;
            unlim_urgent_doses = o.unlim_urgent_doses;

            reactive_vac_strategy = o.reactive_vac_strategy;
            reactive_vac_dose_allocation = o.reactive_vac_dose_allocation;

            min_age = o.min_age;
        }

        virtual ~Vac_Campaign() { // is this correct to address a vtable error?
            for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                while (eligibility_queue[dose].size()) {
                    Eligibility_Group* eg = eligibility_queue[dose].top();
                    eligibility_queue[dose].pop();
                    delete eg;
                }

                while (urgent_queue[dose].size()) {
                    Eligibility_Group* eg = urgent_queue[dose].top();
                    urgent_queue[dose].pop();
                    delete eg;
                }
            }
        }

        void set_par(const Parameters* par) { _par = par; }
        void set_rng(gsl_rng* rng) { _VAX_RNG = rng; }

        int get_start_of_campaign(int campaignType) { return start_of_campaign[campaignType]; }
        void set_start_of_campaign(int campaignType, int day) { start_of_campaign[campaignType] = day; }

        int get_end_of_campaign(int campaignType) { return end_of_campaign[campaignType]; }
        void set_end_of_campaign(int campaignType, int day) { end_of_campaign[campaignType] = day; }

        void set_prioritize_first_doses(bool val)    { prioritize_first_doses = val; }
        void set_flexible_queue_allocation(bool val) { flexible_queue_allocation = val; }
        void set_unlim_urgent_doses(bool val)        { unlim_urgent_doses = val; }

        void set_reactive_vac_strategy(VacCampaignType vct) { reactive_vac_strategy = vct; }
        VacCampaignType get_reactive_vac_strategy() { return reactive_vac_strategy; }

        void set_reactive_vac_dose_allocation(double val) { reactive_vac_dose_allocation = val; }
        double get_reactive_vac_dose_allocation() { return reactive_vac_dose_allocation; }

        std::vector< std::priority_queue<Eligibility_Group*, std::vector<Eligibility_Group*>, Eligibility_Group::comparator> > get_eligibility_queue() const { return eligibility_queue; }
        void set_eligibility_queue(std::vector< std::priority_queue<Eligibility_Group*, std::vector<Eligibility_Group*>, Eligibility_Group::comparator> > eq) { eligibility_queue = eq; }

        std::vector< std::priority_queue<Eligibility_Group*, std::vector<Eligibility_Group*>, Eligibility_Group::comparator> > get_urgent_queue() const { return urgent_queue; }
        void set_urgent_queue(std::vector< std::priority_queue<Eligibility_Group*, std::vector<Eligibility_Group*>, Eligibility_Group::comparator> > uq) { urgent_queue = uq; }

        std::vector< std::vector< std::map<int, int> > > get_orig_doses_available() const { return orig_doses_available; }
        int get_doses_available(int day, int dose, int age_bin) { return doses_available[day][dose][age_bin]; }
        int get_doses_available(int day) const {
            int tot_doses = 0;
            for (auto vec : doses_available[day]) {
                for (auto const& [bin, doses] : vec) {
                    tot_doses += doses;
                }
            }
            return tot_doses;
        }

        void set_urg_doses_available(std::vector< std::vector< std::map<int, int> > > da) { urg_doses_available = da; }
        void set_doses_available(std::vector< std::vector< std::map<int, int> > > da) { doses_available = da; }
        void set_doses_available(int day, int dose, int age_bin, int da) { doses_available[day][dose][age_bin] = da; }

        void init_orig_doses_available() {
            for (int day = 0; day < (int) _par->runLength; ++day) {
                for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                    for (int bin : unique_age_bins) {
                        orig_doses_available[day][dose][bin] = doses_available[day][dose][bin] + urg_doses_available[day][dose][bin];
                    }
                }
            }
        }

        void _adj_orig_doses(int day, int dose, int age_bin, int da) {
            orig_doses_available[day+1][dose][age_bin] += da;
            orig_doses_available[day][dose][age_bin] -= da;
        }

        void rollover_unused_doses(int day, int dose, int age_bin) {
            int leftover_std_doses = doses_available[day][dose][age_bin];
            int leftover_urg_doses = urg_doses_available[day][dose][age_bin];

            if (((day + 1) < (int) _par->runLength) and (leftover_std_doses + leftover_urg_doses)) {
                doses_available[day][dose][age_bin] = 0;
                doses_available[day+1][dose][age_bin] += leftover_std_doses;

                urg_doses_available[day][dose][age_bin] = 0;
                urg_doses_available[day+1][dose][age_bin] += leftover_urg_doses;

                _adj_orig_doses(day, dose, age_bin, leftover_std_doses + leftover_urg_doses);
            }
        }

        std::vector<int> get_min_age() const  { return min_age; }
        int get_min_age(size_t day)    const  { return min_age.at(day); }
        void set_min_age(std::vector<int> ma) { min_age = ma; }

        bool is_age_eligible_on(int age, int day) { return age >= min_age.at(day); }

        std::vector< std::map<int, std::vector<Person*> > > get_potential_vaccinees() const { return potential_vaccinees; }
        std::vector<Person*> get_potential_vaccinees(int dose, int age_bin) { return potential_vaccinees[dose][age_bin]; }
        void set_potential_vaccinees(std::vector< std::map<int, std::vector<Person*> > > pv) { potential_vaccinees = pv; }
        void set_potential_vaccinees(int dose, int age_bin, std::vector<Person*> pv) { potential_vaccinees[dose][age_bin] = pv; }

        std::vector< std::map<int, std::vector<Person*> > > get_urgent_vaccinees() const { return urgent_vaccinees; }
        void set_urgent_vaccinees(std::vector< std::map<int, std::vector<Person*> > > uv) { urgent_vaccinees = uv; }

        void add_potential_vaccinee(int dose, int age_bin, Person* p) { potential_vaccinees[dose][age_bin].push_back(p); }

        void add_new_eligible_people(int today) {
            for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {    // only need to iterate for any dose that would require a next dose
                Eligibility_Group* std_eg = nullptr;
                Eligibility_Group* urg_eg = nullptr;
                if ((eligibility_queue[dose].size() > 0) and (eligibility_queue[dose].top()) and (eligibility_queue[dose].top()->eligibility_day == today)) {
                    std_eg = eligibility_queue[dose].top();
                }

                if ((urgent_queue[dose].size() > 0) and (urgent_queue[dose].top()) and (urgent_queue[dose].top()->eligibility_day == today)) {
                    urg_eg = urgent_queue[dose].top();
                }

                for (int bin : unique_age_bins) {
                    if (std_eg) {
                        std::vector<Person*> std_people_to_add = std_eg->eligible_people[bin];
                        potential_vaccinees[dose][bin].insert(potential_vaccinees[dose][bin].end(), std_people_to_add.begin(), std_people_to_add.end());
                    }

                    if (urg_eg) {
                        std::vector<Person*> urg_people_to_add = urg_eg->eligible_people[bin];
                        urgent_vaccinees[dose][bin].insert(urgent_vaccinees[dose][bin].end(), urg_people_to_add.begin(), urg_people_to_add.end());
                    }

                }

                if (std_eg) {
                    eligibility_queue[dose].pop();
                    delete std_eg;
                }

                if (urg_eg) {
                    urgent_queue[dose].pop();
                    delete urg_eg;
                }
            }
        }

        // swap the value at index with the value at the end
        // warning: v is passed by reference and will be altered
        void _move_to_end(int index, std::vector<Person*>& v) {
            assert(index < (int) v.size());
            Person* tmp = v.back();
            v.back() = v[index];
            v[index] = tmp;
        }

        Vaccinee* next_vaccinee(size_t day, int dose, int age_bin) {
            Person* person = nullptr;
            VaccineAllocationType source_of_dose = NUM_OF_VACCINE_ALLOCATION_TYPES;
            VaccinationQueueType  vaccinee_queue = NUM_OF_VACCINATION_QUEUE_TYPES;

            // TODO: handle first dose prioritization and flexible dose allocation

            int std_doses = doses_available[day][dose][age_bin];
            int urg_doses = urg_doses_available[day][dose][age_bin];

            std::vector<Person*>& urgent_pool   = urgent_vaccinees[dose][age_bin];
            std::vector<Person*>& vaccinee_pool = potential_vaccinees[dose][age_bin];

            if (std_doses + urg_doses) {
                std::vector<Person*>* pool = nullptr;
                if (urg_doses and urgent_pool.size()) {
                    pool = &urgent_pool;
                    source_of_dose = URGENT_ALLOCATION;
                    vaccinee_queue = URGENT_QUEUE;
                } else if (std_doses and vaccinee_pool.size()) {
                    pool = &vaccinee_pool;
                    source_of_dose = STANDARD_ALLOCATION;
                    vaccinee_queue = STANDARD_QUEUE;
                }

                if (pool) {
                    int ran_index = gsl_ran_flat(_VAX_RNG, 0, pool->size());     // gsl selects on a<=x<b
                    _move_to_end(ran_index, *pool);
                    person = pool->back();
                }
            }

            Vaccinee* vaccinee = person ? new Vaccinee(person, vaccinee_queue, source_of_dose) : nullptr;
            return vaccinee;
        }

        void remove_from_pool(int dose, int bin, Person* p) {
            assert((not potential_vaccinees[dose][bin].empty()) and (p == potential_vaccinees[dose][bin].back()));
            potential_vaccinees[dose][bin].pop_back();
        }

        bool vaccinate(Vaccinee* v, int day) { return is_age_eligible_on(v->get_person()->getAge(), day) and v->vaccinate(day); }

        void tally_dose(int day, int dose, int age_bin, Vaccinee* v) {
            if (v->get_dose_source() == STANDARD_ALLOCATION) {
                --doses_available[day][dose][age_bin];
                ++doses_used[day][dose][age_bin];
            } else if (v->get_dose_source() == URGENT_ALLOCATION) {
                --urg_doses_available[day][dose][age_bin];
                ++urg_doses_used[day][dose][age_bin];
            }
        }

        int get_doses_used(int day, int dose, int age_bin) { return doses_used[day][dose][age_bin]; }
        int get_urg_doses_used(int day, int dose, int age_bin) { return urg_doses_used[day][dose][age_bin]; }

        void schedule_revaccinations(vector<Eligibility_Group*> revaccinations) {
            assert(eligibility_queue.size() == revaccinations.size());
            for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                if (revaccinations[dose]) { eligibility_queue[dose].push(revaccinations[dose]); }
            }
        }

        void schedule_urgent_doses(vector<Eligibility_Group*> urgents) {
            assert(eligibility_queue.size() == urgents.size());
            for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                if (urgents[dose]) { eligibility_queue[dose].push(urgents[dose]); }
            }
        }

        void ring_scheduling(int day, vector< set<Person*> > tracedContacts);
        void geographic_scheduling(int day, vector< set<Person*> > targetedPeople, Community* community);
        void location_scheduling(int day, vector< set<Person*> > targetedPeople);

        void reactive_strategy(int day, vector< set<Person*> > targetedPeople, Community* community);

        int get_age_bin(int age) const { return age_bin_lookup.at(age); }
        std::vector<int> get_unique_age_bins() const { return unique_age_bins; }
        std::map<int, int> get_unique_age_bin_pops() const { return unique_age_bin_pops; }

        void generate_age_bins(Community* community, std::set<int> unique_bin_mins, std::set<int> unique_bin_maxs);

        void init_eligibility_queue(const Community* community);

        vector<Eligibility_Group*> init_new_eligible_groups(int day);

    private:
        std::vector< std::vector< std::map<int, int> > > orig_doses_available;      // copy of total pre-sim doses_available for use in model analytics
        std::vector< std::vector< std::map<int, int> > > doses_available;           // daily standard dose availability indexed by [day][dose][age bin]
        std::vector< std::vector< std::map<int, int> > > urg_doses_available;       // daily urgent dose availability indexed by [day][dose][age bin]

        std::vector< std::vector< std::map<int, int> > > doses_used;                // daily doses used indexed by [day][dose][age bin]
        std::vector< std::vector< std::map<int, int> > > urg_doses_used;                // daily doses used indexed by [day][dose][age bin]

        std::vector< std::map<int, std::vector<Person*> > > potential_vaccinees;    // pool of potential people to be vaccinated indexed by [next dose to give][age bin]
        std::vector< std::map<int, std::vector<Person*> > > urgent_vaccinees;       // pool of potential people to be urgently vaccinated indexed by [next dose to give][age bin]

        // outer vector indexed by next dose dose to be given
        // each queue holds the list of eligibility groups ordered by date of next dose
        std::vector< std::priority_queue<Eligibility_Group*, std::vector<Eligibility_Group*>, Eligibility_Group::comparator> > eligibility_queue;
        std::vector< std::priority_queue<Eligibility_Group*, std::vector<Eligibility_Group*>, Eligibility_Group::comparator> > urgent_queue;

        std::vector<int> age_bin_lookup;
        std::vector<int> unique_age_bins;
        std::map<int, int> unique_age_bin_pops;

        const Parameters* _par;
        gsl_rng* _VAX_RNG;

        std::vector<int> start_of_campaign;                    // index by VacCampaignType and returns start date as sim day
        std::vector<int> end_of_campaign;                      // index by VacCampaignType and returns start date as sim day

        bool prioritize_first_doses;                           // do we prioritize first doses, or completing vac schedules? --> possibly move to _par?
        bool flexible_queue_allocation;                        // can doses be used from any allocation, or only as intended?
        bool unlim_urgent_doses;                               // should an unlimited number of doses be used for the urgent queue (i.e. reactive strategy)?

        bool pool_urg_doses;                                   // urgent doses are pooled and used for any urgent vaccinee (regardless of age or dose)
        bool pool_std_doses;                                   // standard doses are pooled and used for any standard vaccinee (regardless of age or dose)
        bool pool_all_doses;                                   // all doses are pooled and used for any vaccinee (regardless of age or dose)

        VacCampaignType reactive_vac_strategy;                 // parameter for type of reactive strategy (if one is active)
        double reactive_vac_dose_allocation;                   // what proportion of total daily doses are reserved for reactive strategies

        std::vector<int> min_age;
};


#endif
