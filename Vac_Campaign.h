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

// will hold people scheduled to become eleigible for a certain vaccine dose
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

// holds pre-constructed Eligibility_Groups by [dose]
typedef std::vector< std::priority_queue<Eligibility_Group*, std::vector<Eligibility_Group*>, Eligibility_Group::comparator> > Eligibility_Q;
// pools of eleigible people to be vaccinated by [dose][age bin]
typedef std::vector< std::map<int, std::vector<Person*> > > Vaccinee_Pool;
// structured storage of pointers to doses available (to allow for pooling) by [day][dose][age bin]
typedef std::vector< std::vector< std::map<int, int*> > > Dose_Ptrs;
// structured storages of values for doses used by [day][dose][age bin]
typedef std::vector< std::vector< std::map<int, int> > > Dose_Vals;

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

            pool_urg_doses = false;
            pool_std_doses = false;
            pool_all_doses = false;

            start_of_campaign = vector<int>(NUM_OF_VAC_CAMPAIGN_TYPES, 0);
            end_of_campaign = vector<int>(NUM_OF_VAC_CAMPAIGN_TYPES, 0);

            std_doses_available = Dose_Ptrs(par->runLength, std::vector< std::map<int, int*> >(par->numVaccineDoses));
            urg_doses_available = Dose_Ptrs(par->runLength, std::vector< std::map<int, int*> >(par->numVaccineDoses));

            std_doses_used = Dose_Vals(par->runLength, std::vector< std::map<int, int> >(par->numVaccineDoses));
            urg_doses_used = Dose_Vals(par->runLength, std::vector< std::map<int, int> >(par->numVaccineDoses));

            potential_std_vaccinees = Vaccinee_Pool(par->numVaccineDoses);
            potential_urg_vaccinees = Vaccinee_Pool(par->numVaccineDoses);

            std_eligibility_queue.clear();
            std_eligibility_queue.resize(_par->numVaccineDoses);

            urg_eligibility_queue.clear();
            urg_eligibility_queue.resize(_par->numVaccineDoses);
        }

        // copy ctor does member-wise copying that will be updated in the Community copy ctor
        Vac_Campaign(const Vac_Campaign& o) {
            _doses = o._doses;
            std_doses_available = o.std_doses_available;
            urg_doses_available = o.urg_doses_available;

            std_doses_used = o.std_doses_used;
            urg_doses_used = o.urg_doses_used;

            potential_std_vaccinees = o.potential_std_vaccinees;
            potential_urg_vaccinees = o.potential_urg_vaccinees;

            std_eligibility_queue = o.std_eligibility_queue;
            urg_eligibility_queue = o.urg_eligibility_queue;

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

            pool_urg_doses = o.pool_urg_doses;
            pool_std_doses = o.pool_std_doses;
            pool_all_doses = o.pool_all_doses;

            reactive_vac_strategy = o.reactive_vac_strategy;
            reactive_vac_dose_allocation = o.reactive_vac_dose_allocation;

            min_age = o.min_age;
        }

        virtual ~Vac_Campaign() { // is this correct to address a vtable error?
            for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                while (std_eligibility_queue[dose].size()) {
                    Eligibility_Group* eg = std_eligibility_queue[dose].top();
                    std_eligibility_queue[dose].pop();
                    delete eg;
                }

                while (urg_eligibility_queue[dose].size()) {
                    Eligibility_Group* eg = urg_eligibility_queue[dose].top();
                    urg_eligibility_queue[dose].pop();
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
        void set_pool_urg_doses(bool val)            { pool_urg_doses = val; }
        void set_pool_std_doses(bool val)            { pool_std_doses = val; }
        void set_pool_all_doses(bool val)            { pool_all_doses = val; }

        void set_reactive_vac_strategy(VacCampaignType vct) { reactive_vac_strategy = vct; }
        VacCampaignType get_reactive_vac_strategy() { return reactive_vac_strategy; }

        void set_reactive_vac_dose_allocation(double val) { reactive_vac_dose_allocation = val; }
        double get_reactive_vac_dose_allocation() { return reactive_vac_dose_allocation; }

        Eligibility_Q get_urg_eligibility_queue() const { return urg_eligibility_queue; }
        void set_urg_eligibility_queue(Eligibility_Q uq) { urg_eligibility_queue = uq; }

        Eligibility_Q get_std_eligibility_queue() const { return std_eligibility_queue; }
        void set_std_eligibility_queue(Eligibility_Q eq) { std_eligibility_queue = eq; }

        void init_doses_available(Dose_Vals urg_in, Dose_Vals std_in);

        int get_urg_doses_available(int day, int dose, int age_bin) { return *urg_doses_available[day][dose][age_bin]; }
        int get_std_doses_available(int day, int dose, int age_bin) { return *std_doses_available[day][dose][age_bin]; }

        Dose_Ptrs get_urg_doses_available() { return urg_doses_available; }
        Dose_Ptrs get_std_doses_available() { return std_doses_available; }

        // aggregates all std and urg doses available for a given day
        int get_all_doses_available(int day) {
            int tot_doses = 0;
            for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                for (int bin : unique_age_bins) {
                    tot_doses += (*std_doses_available[day][dose][bin]) + (*urg_doses_available[day][dose][bin]);
                }
            }
            return tot_doses;
        }

        void set_urg_doses_available(Dose_Ptrs da) { urg_doses_available = da; }
        void set_std_doses_available(Dose_Ptrs da) { std_doses_available = da; }

        void set_urg_doses_available(int day, int dose, int age_bin, int da) { *urg_doses_available[day][dose][age_bin] = da; }
        void set_std_doses_available(int day, int dose, int age_bin, int da) { *std_doses_available[day][dose][age_bin] = da; }

        void set_urg_dose_ptr(int day, int dose, int age_bin, int* ptr) { urg_doses_available[day][dose][age_bin] = ptr; }
        void set_std_dose_ptr(int day, int dose, int age_bin, int* ptr) { std_doses_available[day][dose][age_bin] = ptr; }

        void rollover_unused_doses(int day, int dose, int age_bin) {
            int leftover_std_doses = *std_doses_available[day][dose][age_bin];
            int leftover_urg_doses = *urg_doses_available[day][dose][age_bin];

            if (((day + 1) < (int) _par->runLength) and (leftover_std_doses + leftover_urg_doses)) {
                *std_doses_available[day][dose][age_bin] = 0;
                *std_doses_available[day+1][dose][age_bin] += leftover_std_doses;

                *urg_doses_available[day][dose][age_bin] = 0;
                *urg_doses_available[day+1][dose][age_bin] += leftover_urg_doses;
            }
        }

        std::vector<int> get_min_age() const  { return min_age; }
        int get_min_age(size_t day)    const  { return min_age.at(day); }
        void set_min_age(std::vector<int> ma) { min_age = ma; }

        bool is_age_eligible_on(int age, int day) { return age >= min_age.at(day); }

        Vaccinee_Pool get_potential_std_vaccinees() const { return potential_std_vaccinees; }
        std::vector<Person*> get_potential_std_vaccinees(int dose, int age_bin) { return potential_std_vaccinees[dose][age_bin]; }

        void set_potential_std_vaccinees(Vaccinee_Pool pv) { potential_std_vaccinees = pv; }
        void set_potential_std_vaccinees(int dose, int age_bin, std::vector<Person*> pv) { potential_std_vaccinees[dose][age_bin] = pv; }

        Vaccinee_Pool get_potential_urg_vaccinees() const { return potential_urg_vaccinees; }
        void set_potential_urg_vaccinees(Vaccinee_Pool uv) { potential_urg_vaccinees = uv; }

        void add_potential_std_vaccinee(int dose, int age_bin, Person* p) { potential_std_vaccinees[dose][age_bin].push_back(p); }

        // will be called everyday to check if new people need to be added to the pools
        void add_new_eligible_people(int today) {
            for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                Eligibility_Group* std_eg = nullptr;
                Eligibility_Group* urg_eg = nullptr;

                std_eg = _if_valid_eligibility_group_today(std_eligibility_queue, dose, today);
                urg_eg = _if_valid_eligibility_group_today(urg_eligibility_queue, dose, today);

                for (int bin : unique_age_bins) {
                    if (std_eg) { _insert_eligible_people(std_eg, potential_std_vaccinees, dose, bin); }
                    if (urg_eg) { _insert_eligible_people(urg_eg, potential_urg_vaccinees, dose, bin); }
                }

                // as we sim, Eligibility_Groups will be deleted (some will be left to be cleaned by the dtor)
                if (std_eg) { std_eligibility_queue[dose].pop(); delete std_eg; }
                if (urg_eg) { urg_eligibility_queue[dose].pop(); delete urg_eg; }
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

        // main function that controls who will be vaccinated next
        // implicit urgent > standard heirarchy
        Vaccinee* next_vaccinee(size_t day, int dose, int age_bin) {
            Person* person = nullptr;
            VaccineAllocationType source_of_dose = NUM_OF_VACCINE_ALLOCATION_TYPES;
            VaccinationQueueType  vaccinee_queue = NUM_OF_VACCINATION_QUEUE_TYPES;

            // TODO: handle first dose prioritization and flexible dose allocation

            int std_doses = *std_doses_available[day][dose][age_bin];
            int urg_doses = *urg_doses_available[day][dose][age_bin];

            std::vector<Person*>& urg_pool = potential_urg_vaccinees[dose][age_bin];
            std::vector<Person*>& std_pool = potential_std_vaccinees[dose][age_bin];

            // only continue if there are any doses available for this day, dose, bin combination
            if (std_doses + urg_doses) {
                std::vector<Person*>* pool = nullptr;
                if (urg_doses and urg_pool.size()) {
                    pool = &urg_pool;
                    source_of_dose = URGENT_ALLOCATION;
                    vaccinee_queue = URGENT_QUEUE;
                } else if (std_doses and std_pool.size()) {
                    pool = &std_pool;
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

        void remove_from_pool(int dose, int bin, Vaccinee* v) {
            Person* p = v->get_person();
            if (v->get_dose_source() == STANDARD_ALLOCATION) {
                assert((not potential_std_vaccinees[dose][bin].empty()) and (p == potential_std_vaccinees[dose][bin].back()));
                potential_std_vaccinees[dose][bin].pop_back();
            } else if (v->get_dose_source() == URGENT_ALLOCATION) {
                assert((not potential_urg_vaccinees[dose][bin].empty()) and (p == potential_urg_vaccinees[dose][bin].back()));
                potential_urg_vaccinees[dose][bin].pop_back();
            }
        }

        bool vaccinate(Vaccinee* v, int day) { return is_age_eligible_on(v->get_person()->getAge(), day) and v->vaccinate(day); }

        void tally_dose(int day, int dose, int age_bin, Vaccinee* v) {
            if (v->get_dose_source() == STANDARD_ALLOCATION) {
                --(*std_doses_available[day][dose][age_bin]);
                ++std_doses_used[day][dose][age_bin];
            } else if (v->get_dose_source() == URGENT_ALLOCATION) {
                --(*urg_doses_available[day][dose][age_bin]);
                ++urg_doses_used[day][dose][age_bin];
            }
        }

        Dose_Vals get_std_doses_used() { return std_doses_used; }
        Dose_Vals get_urg_doses_used() { return urg_doses_used; }

        int get_std_doses_used(int day, int dose, int age_bin) { return std_doses_used[day][dose][age_bin]; }
        int get_urg_doses_used(int day, int dose, int age_bin) { return urg_doses_used[day][dose][age_bin]; }

        void schedule_revaccinations(vector<Eligibility_Group*> revaccinations) { _add_new_eligibility_groups(std_eligibility_queue, revaccinations); }
        void schedule_urgent_doses(vector<Eligibility_Group*> urgents)          { _add_new_eligibility_groups(urg_eligibility_queue, urgents); }

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

        void copy_doses_available(Vac_Campaign* vc);

    private:
        std::vector<int> _doses;                                // doses that other data structures will point to (depending on the pooling set by the user)
        Dose_Ptrs std_doses_available;                           // daily standard dose availability indexed by [day][dose][age bin]
        Dose_Ptrs urg_doses_available;                           // daily urgent dose availability indexed by [day][dose][age bin]

        Dose_Vals std_doses_used;                                    // daily doses used indexed by [day][dose][age bin]
        Dose_Vals urg_doses_used;                                // daily doses used indexed by [day][dose][age bin]

        Vaccinee_Pool potential_std_vaccinees;                   // pool of potential people to be vaccinated indexed by [next dose to give][age bin]
        Vaccinee_Pool potential_urg_vaccinees;                   // pool of potential people to be urgently vaccinated indexed by [next dose to give][age bin]

        // outer vector indexed by next dose dose to be given
        // each queue holds the list of eligibility groups ordered by date of next dose
        Eligibility_Q std_eligibility_queue;
        Eligibility_Q urg_eligibility_queue;

        std::vector<int> age_bin_lookup;
        std::vector<int> unique_age_bins;
        std::map<int, int> unique_age_bin_pops;

        const Parameters* _par;
        gsl_rng* _VAX_RNG;

        std::vector<int> start_of_campaign;                     // index by VacCampaignType and returns start date as sim day
        std::vector<int> end_of_campaign;                       // index by VacCampaignType and returns start date as sim day

        bool prioritize_first_doses;                            // do we prioritize first doses, or completing vac schedules? --> possibly move to _par?
        bool flexible_queue_allocation;                         // can doses be used from any allocation, or only as intended?
        bool unlim_urgent_doses;                                // should an unlimited number of doses be used for the urgent queue (i.e. reactive strategy)?

        bool pool_urg_doses;                                    // urgent doses are pooled and used for any urgent vaccinee (regardless of age or dose)
        bool pool_std_doses;                                    // standard doses are pooled and used for any standard vaccinee (regardless of age or dose)
        bool pool_all_doses;                                    // all doses are pooled and used for any vaccinee (regardless of age or dose)

        VacCampaignType reactive_vac_strategy;                  // parameter for type of reactive strategy (if one is active)
        double reactive_vac_dose_allocation;                    // what proportion of total daily doses are reserved for reactive strategies

        std::vector<int> min_age;

        Dose_Ptrs _dose_pool(Dose_Vals doses_in);
        Dose_Ptrs _dose_store(Dose_Vals doses_in);
        // Dose_Ptrs _init_orig_doses();

        int _deref(int* ptr) { return *ptr; }

        void _clear_and_resize_doses();

        Eligibility_Group* _if_valid_eligibility_group_today(Eligibility_Q& eq, const int dose, const int day) {
            return ((eq[dose].size() > 0) and (eq[dose].top()) and (eq[dose].top()->eligibility_day == day)) ? eq[dose].top() : nullptr;
        }

        void _insert_eligible_people(Eligibility_Group*& eg, Vaccinee_Pool& vp, const int dose, const int bin) {
            std::vector<Person*> people_to_add = eg->eligible_people[bin];
            vp[dose][bin].insert(vp[dose][bin].end(), people_to_add.begin(), people_to_add.end());
        }

        void _add_new_eligibility_groups(Eligibility_Q& eq, const vector<Eligibility_Group*>& new_el_groups) {
            assert(eq.size() == new_el_groups.size());
            for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                if (new_el_groups[dose]) { eq[dose].push(new_el_groups[dose]); }
            }
        }
};


#endif
