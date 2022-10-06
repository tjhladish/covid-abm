// Vac_Campaign.h
#ifndef __VAC_CAMPAIGN_H
#define __VAC_CAMPAIGN_H

#include "Person.h"
#include <deque>
#include <map>
#include <queue>
#include <set>
#include <gsl/gsl_randist.h>

class Person;

enum VaccinationQueueType {
    URGENT_QUEUE,
    STANDARD_QUEUE,
    // REVACCINATE_QUEUE,       // outdated
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

enum VacCampaignType { // VacStrategyType
    NO_CAMPAIGN,
    GENERAL_CAMPAIGN,
    RING_VACCINATION,
    GEO_VACCINATION,
    LOCATION_VACCINATION,
    GROUPED_RISK_VACCINATION,
    RISK_VACCINATION,
    NUM_OF_VAC_CAMPAIGN_TYPES
};

inline std::ostream& operator<<(std::ostream& out, const VacCampaignType value){
    const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value){
        PROCESS_VAL(NO_CAMPAIGN);
        PROCESS_VAL(GENERAL_CAMPAIGN);
        PROCESS_VAL(RING_VACCINATION);
        PROCESS_VAL(GEO_VACCINATION);
        PROCESS_VAL(LOCATION_VACCINATION);
        PROCESS_VAL(GROUPED_RISK_VACCINATION);
        PROCESS_VAL(RISK_VACCINATION);
        PROCESS_VAL(NUM_OF_VAC_CAMPAIGN_TYPES);
    }
#undef PROCESS_VAL
    return out << s;
}

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
    std::map<int, std::vector<Person*>> eligible_people;        // people to become eligible indexed by [age bin]

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

//typedef std::set<Person*, PerPtrComp> Person_Set;

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

            prioritize_first_doses       = false;
            flexible_queue_allocation    = false;
            unlim_urgent_doses           = false;
            reactive_vac_strategy        = NO_CAMPAIGN;
            reactive_vac_dose_allocation = 0.0;

            pool_urg_doses               = false;
            pool_std_doses               = false;
            pool_all_doses               = false;

            start_of_campaign = vector<int>(NUM_OF_VAC_CAMPAIGN_TYPES, 0);
            end_of_campaign = vector<int>(NUM_OF_VAC_CAMPAIGN_TYPES, 0);

            doses_available = vector<Dose_Ptrs>(NUM_OF_VACCINE_ALLOCATION_TYPES, Dose_Ptrs(par->runLength, std::vector< std::map<int, int*> >(par->numVaccineDoses)));

            doses_used = vector<Dose_Vals>(NUM_OF_VACCINE_ALLOCATION_TYPES, Dose_Vals(par->runLength, vector< map<int, int> >(par->numVaccineDoses)));

            potential_vaccinees = vector<Vaccinee_Pool>(NUM_OF_VACCINATION_QUEUE_TYPES, Vaccinee_Pool(par->numVaccineDoses));

            eligibility_queue = vector<Eligibility_Q>(NUM_OF_VACCINATION_QUEUE_TYPES);
            for (Eligibility_Q& eq : eligibility_queue) {
                eq.clear();
                eq.resize(_par->numVaccineDoses);
            }
        }

        // copy ctor does member-wise copying that will be updated in the Community copy ctor
        Vac_Campaign(const Vac_Campaign& o) {
            _doses = o._doses;
            doses_available = o.doses_available;

            doses_used = o.doses_used;

            potential_vaccinees = o.potential_vaccinees;

            eligibility_queue = o.eligibility_queue;

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
                for (Eligibility_Q eq : eligibility_queue) {
                    while (eq[dose].size()) {
                        Eligibility_Group* eg = eq[dose].top();
                        eq[dose].pop();
                        delete eg;
                    }
                }
            }
        }

        enum GroupedRiskDef {
            BY_FILE,
            BY_QUANTILE,
            NUM_OF_GROUPED_RISK_DEF_TYPES
        };

        void set_par(const Parameters* par) { _par = par; }
        void set_rng(gsl_rng* rng) { _VAX_RNG = rng; }

        int get_start_of_campaign(int campaignType) { return start_of_campaign[campaignType]; }
        void set_start_of_campaign(int campaignType, int day) { start_of_campaign[campaignType] = day; }

        int get_end_of_campaign(int campaignType) { return end_of_campaign[campaignType]; }
        void set_end_of_campaign(int campaignType, int day) { end_of_campaign[campaignType] = day; }

        bool is_campaign_active(int campaignType, const int day) const { return (day >= start_of_campaign[campaignType]) and (day <= end_of_campaign[campaignType]); }

        bool get_prioritize_first_doses() const   { return prioritize_first_doses; }
        void set_prioritize_first_doses(bool val) { prioritize_first_doses = val; }

        bool get_flexible_queue_allocation() const   { return flexible_queue_allocation; }
        void set_flexible_queue_allocation(bool val) { flexible_queue_allocation = val; }

        bool get_unlim_urgent_doses() const   { return unlim_urgent_doses; }
        void set_unlim_urgent_doses(bool val) { unlim_urgent_doses = val; }

        bool get_pool_urg_doses() const   { return pool_urg_doses; }
        void set_pool_urg_doses(bool val) { pool_urg_doses = val; }

        bool get_pool_std_doses() const   { return pool_std_doses; }
        void set_pool_std_doses(bool val) { pool_std_doses = val; }

        bool get_pool_all_doses() const   { return pool_all_doses; }
        void set_pool_all_doses(bool val) { pool_all_doses = val; }

        GroupedRiskDef get_grouped_risk_def() const   { return grouped_risk_def; }
        void set_grouped_risk_def(GroupedRiskDef val) { grouped_risk_def = val; }

        VacCampaignType get_reactive_vac_strategy() { return reactive_vac_strategy; }
        void set_reactive_vac_strategy(VacCampaignType vct) { reactive_vac_strategy = vct; }

        double get_reactive_vac_dose_allocation() { return reactive_vac_dose_allocation; }
        void set_reactive_vac_dose_allocation(double val) { reactive_vac_dose_allocation = val; }

        Eligibility_Q get_eligibility_queue(VaccinationQueueType vqt) const { return eligibility_queue[vqt]; }
        void set_eligibility_queue(VaccinationQueueType vqt, Eligibility_Q uq) { eligibility_queue[vqt] = uq; }

        int get_risk_quantile_nbins() const { return risk_quantile_nbins; }
        void set_risk_quantile_nbins(int n) { risk_quantile_nbins = n; }

        void init_doses_available(Dose_Vals urg_in, Dose_Vals std_in);

        int get_doses_available(int day, int dose, int age_bin, VaccineAllocationType vat) {
            return *doses_available[vat][day][dose][age_bin];
        }

        Dose_Ptrs get_doses_available(VaccineAllocationType vat) { return doses_available[vat]; }

        // aggregates all std and urg doses available for a given day
        int get_all_doses_available(int day) {
            int tot_doses = 0;

            if (pool_std_doses) {
                tot_doses += *doses_available[STANDARD_ALLOCATION][day][0][unique_age_bins.front()];
            } else {
                for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                    for (int bin : unique_age_bins) {
                        tot_doses += *doses_available[STANDARD_ALLOCATION][day][dose][bin];
                    }
                }
            }

            if (pool_urg_doses) {
                tot_doses += *doses_available[URGENT_ALLOCATION][day][0][unique_age_bins.front()];
            } else {
                for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                    for (int bin : unique_age_bins) {
                        tot_doses += *doses_available[URGENT_ALLOCATION][day][dose][bin];
                    }
                }
            }

            return tot_doses;
        }

        int get_all_doses_available(VaccineAllocationType vat, int day) {
            int tot_doses = 0;

            switch (vat) {
                case URGENT_ALLOCATION: {
                    if (pool_urg_doses) {
                        tot_doses += *doses_available[URGENT_ALLOCATION][day][0][unique_age_bins.front()];
                    } else {
                        for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                            for (int bin : unique_age_bins) {
                                tot_doses += *doses_available[URGENT_ALLOCATION][day][dose][bin];
                            }
                        }
                    }
                    break;
                }
                case STANDARD_ALLOCATION: {
                    if (pool_std_doses) {
                        tot_doses += *doses_available[STANDARD_ALLOCATION][day][0][unique_age_bins.front()];
                    } else {
                        for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                            for (int bin : unique_age_bins) {
                                tot_doses += *doses_available[STANDARD_ALLOCATION][day][dose][bin];
                            }
                        }
                    }
                    break;
                }
                default: { cerr << "ERROR: invalid VaccineAllocationType provided to get_all_doses_available(VaccineAllocationType vat, int day)" << endl; exit(-1); }
            }

            return tot_doses;
        }

        int get_all_doses_available(VaccineAllocationType vat, int day, int dose) {
            int tot_doses = 0;

            switch (vat) {
                case URGENT_ALLOCATION: {
                    if (pool_urg_doses) {
                        tot_doses += *doses_available[URGENT_ALLOCATION][day][dose][unique_age_bins.front()];
                    } else {
                        for (int bin : unique_age_bins) {
                            tot_doses += *doses_available[URGENT_ALLOCATION][day][dose][bin];
                        }
                    }
                    break;
                }
                case STANDARD_ALLOCATION: {
                    if (pool_std_doses) {
                        tot_doses += *doses_available[STANDARD_ALLOCATION][day][dose][unique_age_bins.front()];
                    } else {
                        for (int bin : unique_age_bins) {
                            tot_doses += *doses_available[STANDARD_ALLOCATION][day][dose][bin];
                        }
                    }
                    break;
                }
                default: { cerr << "ERROR: invalid VaccineAllocationType provided to get_all_doses_available(VaccineAllocationType vat, int day, int dose)" << endl; exit(-1); }
            }

            return tot_doses;
        }

        vector< map<int, map<int, int> > > multinomially_distribute_pooled_doses(int day) {
            // if doses are pooled, multinomially sample doses available doses based on the pops of the dose/bin combinations saved in 2D map
            // set up weights for multinomial dose sampling
            const int D = _par->numVaccineDoses;
            const int A = unique_age_bins.size();
            vector< vector<double> > pool_weights(NUM_OF_VACCINATION_QUEUE_TYPES, vector<double>((double)D * A));

            if (get_pool_std_doses()) {
                for (int dose = 0; dose < D; ++dose) {
                    for (int bin = 0; bin < A; ++bin) {
                        Vaccinee_Pool std_pool = get_potential_vaccinees(STANDARD_QUEUE);
                        pool_weights[STANDARD_QUEUE][bin + (dose * A)] = std_pool[dose].count(unique_age_bins[bin]) ? std_pool[dose][unique_age_bins[bin]].size() : 0;
                    }
                }
            }

            if (get_pool_urg_doses()) {
                for (int dose = 0; dose < D; ++dose) {
                    for (int bin = 0; bin < A; ++bin) {
                        Vaccinee_Pool urg_pool = get_potential_vaccinees(URGENT_QUEUE);
                        pool_weights[URGENT_QUEUE][bin + (dose * A)] = urg_pool[dose].count(unique_age_bins[bin]) ? urg_pool[dose][unique_age_bins[bin]].size() : 0;
                    }
                }
            }

            // set up return structure and sample is necessary
            vector< vector<unsigned int> > sampled_doses(NUM_OF_VACCINE_ALLOCATION_TYPES, vector<unsigned int>(D * A));
            vector< map<int, map<int, int> > > daily_sampled_doses_available(NUM_OF_VACCINE_ALLOCATION_TYPES);
            for (int dose = 0; dose < D; ++dose) {
                for (int bin = 0; bin < A; ++bin) {
                    daily_sampled_doses_available[STANDARD_ALLOCATION][dose][unique_age_bins[bin]] = 0;
                    daily_sampled_doses_available[URGENT_ALLOCATION][dose][unique_age_bins[bin]]   = 0;
                }
            }

            if (get_pool_std_doses()) {
                const int N = get_all_doses_available(STANDARD_ALLOCATION, day);
                gsl_ran_multinomial(_VAX_RNG, D * A, N, pool_weights[STANDARD_QUEUE].data(), sampled_doses[STANDARD_ALLOCATION].data());
                for (int dose = 0; dose < D; ++dose) {
                    for (int bin = 0; bin < A; ++bin) {
                        daily_sampled_doses_available[STANDARD_ALLOCATION][dose][unique_age_bins[bin]] = sampled_doses[STANDARD_ALLOCATION][bin + (dose * A)];
                    }
                }
            }

            if (get_pool_urg_doses()) {
                const int N = get_all_doses_available(URGENT_ALLOCATION, day);
                gsl_ran_multinomial(_VAX_RNG, D * A, N, pool_weights[URGENT_QUEUE].data(), sampled_doses[URGENT_ALLOCATION].data());
                for (int dose = 0; dose < D; ++dose) {
                    for (int bin = 0; bin < A; ++bin) {
                        daily_sampled_doses_available[URGENT_ALLOCATION][dose][unique_age_bins[bin]] = sampled_doses[URGENT_ALLOCATION][bin + (dose * A)];
                    }
                }
            }

            return daily_sampled_doses_available;
        }

        void set_doses_available(Dose_Ptrs da, VaccineAllocationType vat) { doses_available[vat] = da; }

        void set_doses_available(int day, int dose, int age_bin, int da, VaccineAllocationType vat) { *doses_available[vat][day][dose][age_bin] = da; }

        void set_dose_ptr(int day, int dose, int age_bin, int* ptr, VaccineAllocationType vat) { doses_available[vat][day][dose][age_bin] = ptr; }

        void rollover_unused_doses(int day, int dose, int age_bin) {
            if ((day + 1) < (int) _par->runLength) {
                const int leftover_std_doses = *doses_available[STANDARD_ALLOCATION][day][dose][age_bin];
                if (leftover_std_doses) {
                    *doses_available[STANDARD_ALLOCATION][day][dose][age_bin] = 0;
                    *doses_available[STANDARD_ALLOCATION][day+1][dose][age_bin] += leftover_std_doses;
                }

                const int leftover_urg_doses = *doses_available[URGENT_ALLOCATION][day][dose][age_bin];
                if (leftover_urg_doses and (day >= start_of_campaign[reactive_vac_strategy])) {
                    *doses_available[URGENT_ALLOCATION][day][dose][age_bin] = 0;
                    *doses_available[URGENT_ALLOCATION][day+1][dose][age_bin] += leftover_urg_doses;
                }
            }
        }

        std::vector<int> get_min_age() const  { return min_age; }
        int get_min_age(size_t day)    const  { return min_age.at(day); }
        void set_min_age(std::vector<int> ma) { min_age = ma; }

        bool is_age_eligible_on(int age, int day) { return age >= min_age.at(day); }

        Vaccinee_Pool get_potential_vaccinees(VaccinationQueueType vqt) const { return potential_vaccinees[vqt]; }
        std::vector<Person*> get_potential_vaccinees(VaccinationQueueType vqt, int dose, int age_bin) const { return potential_vaccinees[vqt][dose].at(age_bin); }

        void set_potential_vaccinees(VaccinationQueueType vqt, Vaccinee_Pool pv) { potential_vaccinees[vqt] = pv; }
        void set_potential_vaccinees(VaccinationQueueType vqt, int dose, int age_bin, std::vector<Person*> pv ) { potential_vaccinees[vqt][dose][age_bin] = pv; }

        void add_potential_vaccinee(VaccinationQueueType vqt, int dose, int age_bin, Person* p) { potential_vaccinees[vqt][dose][age_bin].push_back(p); }

        int get_pool_size(Vaccinee_Pool vp) {
            int tot = 0;
            for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                for (int bin : unique_age_bins) {
                    tot += vp[dose][bin].size();
                }
            }
            return tot;
        }

        int get_pool_size_by_dose(Vaccinee_Pool vp, int dose) {
            int tot = 0;
            for (int bin : unique_age_bins) {
                tot += vp[dose][bin].size();
            }
            return tot;
        }

        // will be called everyday to check if new people need to be added to the pools
        bool add_new_eligible_people(int today) {
            bool group_added = false;
            for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                Eligibility_Group* std_eg = nullptr;
                Eligibility_Group* urg_eg = nullptr;

                std_eg = _if_valid_eligibility_group_today(eligibility_queue[STANDARD_QUEUE], dose, today);
                urg_eg = _if_valid_eligibility_group_today(eligibility_queue[URGENT_QUEUE], dose, today);

                for (int bin : unique_age_bins) {
                    if (std_eg) { _insert_eligible_people(std_eg, potential_vaccinees[STANDARD_QUEUE], dose, bin); }
                    if (urg_eg) { _insert_eligible_people(urg_eg, potential_vaccinees[URGENT_QUEUE], dose, bin); }
                }

                // as we sim, Eligibility_Groups will be deleted (some will be left to be cleaned by the dtor)
                if (std_eg) { eligibility_queue[STANDARD_QUEUE][dose].pop(); delete std_eg; }
                if (urg_eg) { eligibility_queue[URGENT_QUEUE][dose].pop(); delete urg_eg; }

                if (std_eg or urg_eg) { group_added = true; }
            }
            if (reactive_vac_strategy == GROUPED_RISK_VACCINATION and (today > start_of_campaign[GROUPED_RISK_VACCINATION])) {
                _add_ppl_from_risk_groups(today);
            }
            return group_added;
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
        // implicit standard > urgent heirarchy
        Vaccinee* next_vaccinee(size_t day, int dose, int age_bin, const vector< map<int, map<int, int> > >& daily_sampled_doses_available) {
            Person* person = nullptr;
            VaccineAllocationType source_of_dose = NUM_OF_VACCINE_ALLOCATION_TYPES;
            VaccinationQueueType  vaccinee_queue = NUM_OF_VACCINATION_QUEUE_TYPES;

            // TODO: handle first dose prioritization and flexible dose allocation

            const int std_doses = pool_std_doses ? daily_sampled_doses_available[STANDARD_ALLOCATION].at(dose).at(age_bin) : *doses_available[STANDARD_ALLOCATION][day][dose][age_bin];
            const int urg_doses = pool_urg_doses ? daily_sampled_doses_available[URGENT_ALLOCATION].at(dose).at(age_bin) : *doses_available[URGENT_ALLOCATION][day][dose][age_bin];

            std::vector<Person*>& std_pool = potential_vaccinees[STANDARD_QUEUE][dose][age_bin];
            std::vector<Person*>& urg_pool = potential_vaccinees[URGENT_QUEUE][dose][age_bin];

            // only continue if there are any doses available for this day, dose, bin combination
            if (std_doses or urg_doses) {
                std::vector<Person*>* pool = nullptr;
                if (std_doses and std_pool.size()) {
                    pool = &std_pool;
                    source_of_dose = STANDARD_ALLOCATION;
                    vaccinee_queue = STANDARD_QUEUE;
                } else if (urg_doses and urg_pool.size()) {
                    pool = &urg_pool;
                    source_of_dose = URGENT_ALLOCATION;
                    vaccinee_queue = URGENT_QUEUE;
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
            VaccinationQueueType status = v->get_status();

            assert((not potential_vaccinees[status][dose][bin].empty()) and (p == potential_vaccinees[status][dose][bin].back()));
            potential_vaccinees[status][dose][bin].pop_back();
        }

        bool vaccinate(Vaccinee* v, int day) { return is_age_eligible_on(v->get_person()->getAge(), day) and v->vaccinate(day); }

        void tally_dose(int day, int dose, int age_bin, Vaccinee* v, vector< map<int, map<int, int> > >& daily_sampled_doses_available) {
            VaccineAllocationType vat = v->get_dose_source();
            --(*doses_available[vat][day][dose][age_bin]);
            ++doses_used[vat][day][dose][age_bin];

            if (get_pool_std_doses() and (vat == STANDARD_ALLOCATION)) { --daily_sampled_doses_available[STANDARD_ALLOCATION][dose][age_bin]; }
            if (get_pool_urg_doses() and (vat == URGENT_ALLOCATION))   { --daily_sampled_doses_available[URGENT_ALLOCATION][dose][age_bin]; }
        }

        vector<Dose_Vals> get_doses_used() { return doses_used; }
        Dose_Vals get_doses_used(VaccineAllocationType vat) { return doses_used[vat]; }

        int get_doses_used(int day, int dose, int age_bin, VaccineAllocationType vat) { return doses_used[vat][day][dose][age_bin]; }
        int get_doses_used(int day, VaccineAllocationType vat) {
            int tot = 0;
            for (int dose = 0; dose < _par->numVaccineDoses; ++dose) { // TODO - use foreach-style loops over the structure
                for (int bin : unique_age_bins) {
                    tot += doses_used[vat][day][dose][bin];
                }
            }
            return tot;
        }


        int get_all_doses_used(int day, int dose, int age_bin) { return doses_used[URGENT_ALLOCATION][day][dose][age_bin] + doses_used[STANDARD_ALLOCATION][day][dose][age_bin]; }
        int get_all_doses_used(int day) {
            int tot = 0;
            for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                for (int bin : unique_age_bins) {
                    tot += get_all_doses_used(day, dose, bin);
                }
            }
            return tot;
        }

        int get_doses_used_by_dose(int day, int dose, VaccineAllocationType vat) {
            int tot = 0;
            for (int bin : unique_age_bins) {
                tot += doses_used[vat][day][dose][bin];
            }
            return tot;
        }

        void assign_vaccinee_for_revaccination(Vaccinee* v, int dose, int bin, vector<Eligibility_Group*> &std_revac, vector<Eligibility_Group*> &urg_revac) {
            VaccinationQueueType q = v->get_status();
            if (q == URGENT_QUEUE) {
                urg_revac[dose + 1]->eligible_people[bin].push_back(v->get_person());
            } else if (q == STANDARD_QUEUE) {
                std_revac[dose + 1]->eligible_people[bin].push_back(v->get_person());
            } else {
                cerr << "Vacciee not from any queue." << endl;
                exit(-1);
            }
        }

        void schedule_revaccinations(vector<Eligibility_Group*> urg_revaccinations, vector<Eligibility_Group*> std_revaccinations) {
            _add_new_eligibility_groups(eligibility_queue[URGENT_QUEUE], urg_revaccinations);
            _add_new_eligibility_groups(eligibility_queue[STANDARD_QUEUE], std_revaccinations);
        }
        void schedule_urgent_doses(vector<Eligibility_Group*> urgents)          { _add_new_eligibility_groups(eligibility_queue[URGENT_QUEUE], urgents); }

        void ring_scheduling(int day, vector<set<Person*, PerPtrComp>> tracedContacts);
        void geographic_scheduling(int day, vector<set<Person*, PerPtrComp>> targetedPeople, Community* community);
        void location_scheduling(int day, vector<set<Person*, PerPtrComp>> targetedPeople);
        void grouped_risk_scheduling(int day, Community* community);
        void risk_scheduling(int day, Community* community, double hosp_risk_threshold = 0.1); // hrt is Pr{hospitalization|infection} @ pandemic start

        void reactive_strategy(int day, vector<set<Person*, PerPtrComp>> targetedPeople, Community* community);

        int get_age_bin(int age) const { return age_bin_lookup.at(age); }
        std::vector<int> get_unique_age_bins() const { return unique_age_bins; }
        std::map<int, int> get_unique_age_bin_pops() const { return unique_age_bin_pops; }

        void generate_age_bins(Community* community, std::set<int> unique_bin_mins, std::set<int> unique_bin_maxs);

        void init_eligibility_queue(const Community* community);

        vector<Eligibility_Group*> init_new_eligible_groups(int day);

        void copy_doses_available(Vac_Campaign* vc);

        // can be used to create a chace and read from cache
        Vac_Campaign* quick_cache();

        // lookup whether a certain strategy needs contact tracing
        bool contact_tracing_required(VacCampaignType vct) {
            switch (vct) {
                case RING_VACCINATION:          [[fallthrough]];
                case GEO_VACCINATION:           [[fallthrough]];
                case LOCATION_VACCINATION:      return true;
                case NO_CAMPAIGN:               [[fallthrough]];
                case GENERAL_CAMPAIGN:          [[fallthrough]];
                case GROUPED_RISK_VACCINATION:  [[fallthrough]];
                case RISK_VACCINATION:          return false;
                case NUM_OF_VAC_CAMPAIGN_TYPES: [[fallthrough]];
                default:                        cerr << "not valid VacCampaignType" << endl; exit(-1);
            }
        }

    private:
        std::vector<int> _doses;                        // doses that other data structures will point to (depending on the pooling set by the user)
        vector<Dose_Ptrs> doses_available;              // daily dose availability indexed by [VaccineAllocationType][day][dose][age bin]
        vector<Dose_Vals> doses_used;                   // daily doses used indexed by [allocation][day][dose][age bin]

        vector<Vaccinee_Pool> potential_vaccinees;      // pool of potential people to be vaccinated indexed by [VaccinationQueueType][next dose to give][age bin]

        // outer vector indexed by [VaccinationQueueType][next dose dose to be given]
        // each queue holds the list of eligibility groups ordered by date of next dose
        vector<Eligibility_Q> eligibility_queue;

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

        GroupedRiskDef grouped_risk_def;                        // will risk groups be provided by file or generated from quantiles

        VacCampaignType reactive_vac_strategy;                  // parameter for type of reactive strategy (if one is active)
        double reactive_vac_dose_allocation;                    // what proportion of total daily doses are reserved for reactive strategies

        std::vector<int> min_age;

        Dose_Ptrs _dose_pool(Dose_Vals doses_in);
        Dose_Ptrs _dose_store(Dose_Vals doses_in);
        void _dose_pool_all(vector<Dose_Ptrs> &doses_out, Dose_Vals urg_in, Dose_Vals std_in);
        // Dose_Ptrs _init_orig_doses();

        // specialty members for risk group vax strategy
        std::deque< pair<int, vector<Eligibility_Group*>> > _grouped_risk_deque;    // order of groups to be added to the pool
        int _current_risk_group;                                                    // which group was last added to the pool
        std::map<int, vector<Person*>> _sch_risk_groups;                            // people in each group

        int _deref(int* ptr) { return *ptr; }

        void _clear_and_resize_doses();

        Eligibility_Group* _if_valid_eligibility_group_today(Eligibility_Q& eq, const int dose, const int day) {
            return ((eq[dose].size() > 0) and (eq[dose].top()) and (eq[dose].top()->eligibility_day <= day)) ? eq[dose].top() : nullptr;
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

        // specialty method for group risk strategy
        // handles evaluating the stopping criteria for when to move from one group to the next
        bool _ready_to_add_next_group(VaccinationQueueType vqt, VaccineAllocationType vat, int today) {
            return get_pool_size_by_dose(potential_vaccinees[vqt], 0) < (get_all_doses_available(vat, today)); //<= (_sch_risk_groups[_current_risk_group].size() * 0.1);
        }

        // specialty method for group risk strategy
        // adds the next group from the deque if the stopping criteria was met
        void _add_ppl_from_risk_groups(int today) {
            if (_ready_to_add_next_group(URGENT_QUEUE, URGENT_ALLOCATION, today) and _grouped_risk_deque.size()) {
                for (int dose = 0; dose < _par->numVaccineDoses; ++dose) {
                    for (int bin : unique_age_bins) {
                        _insert_eligible_people(_grouped_risk_deque.front().second[dose], potential_vaccinees[URGENT_QUEUE], dose, bin);
                    }
                    delete _grouped_risk_deque.front().second[dose];
                }
                _current_risk_group =  _grouped_risk_deque.front().first;
                _grouped_risk_deque.pop_front();
            }
        }

        void generate_risk_quantiles(Community* community, map<int, vector<Person*>>& grouped_ppl, map<int, double>& grouped_risk, size_t nbin);

        int risk_quantile_nbins;
};


#endif
