// Vac_Campaign.h
#ifndef __VAC_CAMPAIGN_H
#define __VAC_CAMPAIGN_H

#include "Person.h"
#include <deque>
#include <map>
#include <queue>

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

// REFACTOR
class Eligibility_Group {
  public:
    int eligibility_day;                                    // when these people become eligible for vaccination
    std::map<int, std::vector<Person*>> eligible_people;    // people to become eligible index by [age bin]

    struct comparator { bool operator()(const Eligibility_Group* A, const Eligibility_Group* B) const { return A->eligibility_day < B->eligibility_day; } };
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
        Vac_Campaign() {
            prioritize_first_doses = false;
            flexible_queue_allocation = false;
            unlim_urgent_doses = false;
            reactive_vac_strategy = NUM_OF_VAC_CAMPAIGN_TYPES;
            reactive_vac_dose_allocation = 0.0;
            start_of_campaign = vector<int>(NUM_OF_VAC_CAMPAIGN_TYPES, 0);
            end_of_campaign = vector<int>(NUM_OF_VAC_CAMPAIGN_TYPES, 0);
        }
        Vac_Campaign(const Vac_Campaign& o) {
            // urgent_queue                 = o.urgent_queue;
            // standard_queue               = o.standard_queue;
            // revaccinate_queue            = o.revaccinate_queue;
            prioritize_first_doses       = o.prioritize_first_doses;
            doses_available              = o.doses_available;
            // queue_tally                  = o.queue_tally;
            flexible_queue_allocation    = o.flexible_queue_allocation;
            reactive_vac_strategy        = o.reactive_vac_strategy;
            reactive_vac_dose_allocation = o.reactive_vac_dose_allocation;
            min_age                      = o.min_age;
        }
        virtual ~Vac_Campaign() {} //is this correct to address a vtable error?

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

        // void set_revaccination_follow_up(double val) { revaccination_follow_up = val; }

        // size_t get_urgent_queue_size() { return urgent_queue.size(); }
        // size_t get_standard_queue_size() { return standard_queue.size(); }
        // size_t get_revaccinate_queue_size(int day) { return revaccinate_queue[day].size(); }

        // void prioritize_vaccination(Person* p) { urgent_queue.push_back(p); }
        // void schedule_vaccination(Person* p) { standard_queue.push_back(p); }

        // std::deque<Person*> getUrgentQueue() const { return urgent_queue; }
        // void setUrgentQueue(std::deque<Person*> uq) { urgent_queue = uq; }

        // std::deque<Person*> getStandardQueue() const { return standard_queue; }
        // void setStandardQueue(std::deque<Person*> sq) { standard_queue = sq; }

        // std::vector< std::set<Person*, Person::PerPtrComp> > getRevaccinateQueue() const { return revaccinate_queue; }
        // void setRevaccinateQueue(std::vector< std::set<Person*, Person::PerPtrComp> > rq) { revaccinate_queue = rq; }

        // int get_doses_available(int day, VaccineAllocationType alloc) { return doses_available[day][alloc]; }
        // void set_doses_available(std::vector< std::vector<int> > da) {
        //     doses_available = da;
        //     //doses_used.resize(da.size(), std::vector<int>(NUM_OF_VACCINATION_QUEUE_TYPES));
        //     revaccinate_queue.resize(da.size());
        //     queue_tally.resize(da.size(), vector<int>(NUM_OF_VACCINATION_QUEUE_TYPES));
        // }
        int get_doses_available(int day, int dose, int age_bin) { return doses_available[day][dose-1][age_bin]; }
        int get_doses_available(int day) const {
            int tot_doses = 0;
            for (auto vec : doses_available[day]) {
                for (auto const& [bin, doses] : vec) {
                    tot_doses += doses;
                }
            }
            return tot_doses;
        }

        void set_doses_available(std::vector< std::vector< std::map<int, int> > > da) { doses_available = da; }
        void set_doses_available(int day, int dose, int age_bin, int da) { doses_available[day][dose-1][age_bin] = da; }

        std::vector<int> get_min_age() const  { return min_age; }
        int get_min_age(size_t day)    const  { return min_age.at(day); }
        void set_min_age(std::vector<int> ma) { min_age = ma; }

        bool is_age_eligible_on(int age, int day) { return age >= min_age.at(day); }

        // std::map<int, std::vector< std::vector<Person*> >> potential_vaccinees; // pool of potential people to be vaccinated indexed by [age bin min][dose]
        std::vector< std::map<int, std::vector<Person*> > > get_potential_vaccinees() const { return potential_vaccinees; }
        std::vector<Person*> get_potential_vaccinees(int dose, int age_bin) { return potential_vaccinees[dose-1][age_bin]; }

        void set_potential_vaccinees(std::vector< std::map<int, std::vector<Person*> > > pv) { potential_vaccinees = pv; }
        void set_potential_vaccinees(int dose, int age_bin, std::vector<Person*> pv) { potential_vaccinees[dose-1][age_bin] = pv; }
        void add_potential_vaccinee(int dose, int age_bin, Person* pv) { potential_vaccinees[dose-1][age_bin].push_back(pv); }

        void add_new_eligible_people(int today) {
            for (int dose = 0; dose < (_par->numVaccineDoses - 1); ++dose) {    // only need to iterate for any dose that would require a next dose
                if (eligibility_queue[dose].top()->eligibility_day == today) {
                    Eligibility_Group* eg = eligibility_queue[dose].top();
                    for (int bin : unique_age_bins) {  // need to properly implement
                            std::vector<Person*> people_to_add = eg->eligible_people[bin];
                            potential_vaccinees[dose + 1][bin].insert(potential_vaccinees[dose + 1][bin].end(), people_to_add.begin(), people_to_add.end());
                    }
                    eligibility_queue[dose].pop();
                    delete eg;
                }
            }
        }
        /* REFACTOR next_vaccinee()
            check if any doses remain today
            iterate through doses_available to find a dose
            randomly select vaccinee that matches dose reqs from potential_vaccinees
        */
        Vaccinee* next_vaccinee(size_t day, int dose, int age_bin) {
            Person* person = nullptr;
            VaccineAllocationType source_of_dose = NUM_OF_VACCINE_ALLOCATION_TYPES;
            VaccinationQueueType  vaccinee_queue = NUM_OF_VACCINATION_QUEUE_TYPES;

            // TODO: handle urgent/active vaccination here
            // TODO: handle first dose prioritization and flexible dose allocation

            int num_doses_available = doses_available[day][dose][age_bin];
            std::vector<Person*> vaccinee_pool = potential_vaccinees[dose][age_bin];

            if (num_doses_available and vaccinee_pool.size()) {
                int ran_index = gsl_ran_flat(_VAX_RNG, 0, vaccinee_pool.size());     // gsl selects on a<=x<b
                person = vaccinee_pool[ran_index];
                source_of_dose = STANDARD_ALLOCATION;
                vaccinee_queue = STANDARD_QUEUE;
            }

            Vaccinee* vaccinee = person ? new Vaccinee(person, vaccinee_queue, source_of_dose) : nullptr;
            return vaccinee;
        }

        // Vaccinee* next_vaccinee(size_t day) {
        //     // Set up defaults
        //     Person* person = nullptr;
        //     VaccineAllocationType source_of_dose = NUM_OF_VACCINE_ALLOCATION_TYPES;
        //     VaccinationQueueType  vaccinee_queue = NUM_OF_VACCINATION_QUEUE_TYPES;
        //
        //     // Do we have doses available today?
        //     const size_t urgent_doses_avail      = unlim_urgent_doses ? 1 : doses_available.at(day)[URGENT_ALLOCATION];
        //     const size_t standard_doses_avail    = doses_available.at(day)[STANDARD_ALLOCATION];
        //
        //     if (urgent_doses_avail + standard_doses_avail == 0) { return nullptr; } // might as well bail
        //
        //     // How many in each queue are waiting to be vaccinated today?
        //     const size_t urgent_size = urgent_queue.size();
        //     const size_t std_size = standard_queue.size();
        //     const size_t revac_size = revaccinate_queue[day].size();
        //
        //     // TODO -- we may want to have different priorities for finishing a multi-dose vaccine course versus
        //     // administering (e.g. annual) booster doses
        //     // Possible TODO -- this does not currently handle the scenario where the priority is urgent, then
        //     // revaccinations, then standard vaccinations.
        //     if (revac_size and (not prioritize_first_doses or not (urgent_size or std_size))) {
        //     // if (revac_size and (not prioritize_first_doses and not (urgent_size or std_size))) {
        //         // we have someone to revaccinate, and either we're prioritizing completing courses
        //         // or we have no one neededing a first dose
        //         source_of_dose = standard_doses_avail ? STANDARD_ALLOCATION
        //                          : flexible_queue_allocation and urgent_doses_avail ? URGENT_ALLOCATION
        //                          : NUM_OF_VACCINE_ALLOCATION_TYPES;
        //
        //         if (source_of_dose < NUM_OF_VACCINE_ALLOCATION_TYPES) { // i.e., urgent or standard
        //             person = *(revaccinate_queue[day].begin());
        //             revaccinate_queue[day].erase(person);
        //             vaccinee_queue = REVACCINATE_QUEUE;
        //         }
        //     }
        //
        //     if (urgent_size and not person) {
        //         source_of_dose = urgent_doses_avail ? URGENT_ALLOCATION
        //                          : flexible_queue_allocation and standard_doses_avail ? STANDARD_ALLOCATION
        //                          : NUM_OF_VACCINE_ALLOCATION_TYPES;
        //
        //         if (source_of_dose < NUM_OF_VACCINE_ALLOCATION_TYPES) {
        //             person = urgent_queue.front();
        //             urgent_queue.pop_front();
        //             vaccinee_queue = URGENT_QUEUE;
        //         }
        //         if(unlim_urgent_doses) { doses_available.at(day)[URGENT_ALLOCATION]++; }
        //     }
        //
        //     if (std_size and not person) {
        //         source_of_dose = standard_doses_avail ? STANDARD_ALLOCATION
        //                          : flexible_queue_allocation and urgent_doses_avail ? URGENT_ALLOCATION
        //                          : NUM_OF_VACCINE_ALLOCATION_TYPES;
        //
        //         if (source_of_dose < NUM_OF_VACCINE_ALLOCATION_TYPES) {
        //             person = standard_queue.front();
        //             standard_queue.pop_front();
        //             vaccinee_queue = STANDARD_QUEUE;
        //         }
        //     }
        //
        //     // these lines could be combined, but I think the intent is clearer this way
        //     Vaccinee* vaccinee = person ? new Vaccinee(person, vaccinee_queue, source_of_dose) : nullptr;
        //     return vaccinee;
        // }

        bool vaccinate(Vaccinee* v, int day) { return is_age_eligible_on(v->get_person()->getAge(), day) and v->vaccinate(day); }

        // void tally_dose(size_t day, Vaccinee* vac) {
        //     doses_available.at(day)[vac->dose_source]--;
        //     queue_tally[day][vac->status]++;
        // }
        // int get_dose_tally(size_t day, VaccinationQueueType q) { return queue_tally[day][q]; }
        void tally_dose(int day, Vaccinee* vac) {
            int age_bin = age_bin_lookup[vac->get_person()->getAge()];
            int dose = vac->get_person()->getNumVaccinations();

            --doses_available[day][dose-1][age_bin];
            ++doses_used[day][dose-1][age_bin];
        }
        int get_doses_used(int day, int dose, int age_bin) { return doses_used[day][dose-1][age_bin]; }

        void schedule_revaccinations(vector<Eligibility_Group*> revaccinations) {
            for (int dose = 0; dose < (_par->numVaccineDoses - 1); ++dose) {
                eligibility_queue[dose].push(revaccinations[dose]);
            }
        }

        // void schedule_revaccination(size_t day, Vaccinee* vac) {
        //     if (revaccinate_queue.size() > day) {
        //         revaccinate_queue[day].insert(vac->person);
        //         if(vac->get_status() == URGENT_QUEUE and unlim_urgent_doses) {
        //             doses_available.at(day)[STANDARD_ALLOCATION]++;
        //         }
        //     }
        // }

        // void reschedule_remaining_revaccinations(size_t today) {
        //     const size_t queue_length_in_days = revaccinate_queue.size(); // how many days long is the simulation
        //     std::set<Person*, Person::PerPtrComp> remaining_people = revaccinate_queue[today];
        //     if (remaining_people.size() and queue_length_in_days > today + 1) {
        //         revaccinate_queue[today + 1].insert(remaining_people.begin(), remaining_people.end());
        //         revaccinate_queue[today].clear();
        //     }
        // }

        void ring_scheduling(int day, vector< set<Person*> > tracedContacts);
        void geographic_scheduling(int day, vector< set<Person*> > targetedPeople, Community* community);
        void location_scheduling(int day, vector< set<Person*> > targetedPeople);

        void reactive_strategy(int day, vector< set<Person*> > targetedPeople, Community* community);

        int get_age_bin(int age) const { return age_bin_lookup.at(age); }
        std::vector<int> get_unique_age_bins() const { return unique_age_bins; }

        void generate_age_bins(std::set<int> unique_bin_mins, std::set<int> unique_bin_maxs) {
            std::vector<int> mins(unique_bin_mins.begin(), unique_bin_mins.end());
            std::vector<int> maxs(unique_bin_maxs.begin(), unique_bin_maxs.end());

            std::set<int> ages_not_binned;
            assert(mins.size() == maxs.size());

            for (int age = 0; age < NUM_AGE_CLASSES; ++age) {
                ages_not_binned.insert(age);
                for (int i = 0; i < (int)mins.size(); ++i) {
                    if ((age >= mins[i]) and (age <= maxs[i])) {
                        ages_not_binned.erase(age);
                        break;
                    }
                }
            }

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
            age_bin_lookup = std::vector<int>(NUM_AGE_CLASSES);
            for (int age = 0; age < NUM_AGE_CLASSES; ++age) {
                for (int i = 0; i < (int)mins.size(); ++i) {
                    if ((age >= mins[i]) and (age <= maxs[i])) {
                        age_bin_lookup[age] = mins[i];
                        break;
                    }
                }
            }
            unique_age_bins = mins;
        }

        void init_eligibility_queue(const Community* community);

        /*
        void get_reactive_vac_strategy() {
            // check what strategy is being used
            // execute appropriate algorithms to detect and add necessary people to urgent queue

            if(reactive_vac_strategy == NUM_OF_REACTIVE_VAC_STRATEGY_TYPES) { return; }

            switch(reactive_vac_strategy) {
                case RING_VACCINATION:
                    _ring_vaccination();
                    break;
                case GEOGRAPHIC_CLUSTER_VACCINATION:
                    _geographic_cluster_vaccination();
                    break;
            }
        }
        */

    private:
        // REFACTOR
        std::vector< std::vector< std::map<int, int> > > doses_available;         // daily dose availability indexed by [day][dose][age bin]
        std::vector< std::vector< std::map<int, int> > > doses_used;              // daily doses used indexed by [day][dose][age bin]

        std::vector< std::map<int, std::vector<Person*> > > potential_vaccinees; // pool of potential people to be vaccinated indexed by [next dose to give][age bin]

        // outer vector indexed by last dose given
        // each queue holds the list of eligibility groups ordered by date of next dose
        std::vector< std::priority_queue<Eligibility_Group*, std::vector<Eligibility_Group*>, Eligibility_Group::comparator> > eligibility_queue; // TODO: to be tested

        std::vector<int> age_bin_lookup;
        std::vector<int> unique_age_bins;

        const Parameters* _par;
        gsl_rng* _VAX_RNG;

        // std::deque<Person*> urgent_queue;                      // people who need to be vaccinated, determined during simulation
        // std::deque<Person*> standard_queue;                    // people who will be vaccinated, known before transmission sim begins
        // std::vector< std::set<Person*, Person::PerPtrComp> > revaccinate_queue;    // indexed by sim day

        std::vector<int> start_of_campaign;                    // index by VacCampaignType and returns start date as sim day
        std::vector<int> end_of_campaign;                      // index by VacCampaignType and returns start date as sim day

        bool prioritize_first_doses;                           // do we prioritize first doses, or completing vac schedules? --> possibly move to _par?
        bool flexible_queue_allocation;                        // can doses be used from any allocation, or only as intended?
        bool unlim_urgent_doses;                               // should an unlimited number of doses be used for the urgent queue (i.e. reactive strategy)?

        // std::vector< std::vector<int> > doses_available;       // for each day, number of doses available for each queue
        // std::vector< std::vector<int> > queue_tally;           // for each day, number of doses used for each queue

        VacCampaignType reactive_vac_strategy;                 // parameter for type of reactive strategy (if one is active)
        double reactive_vac_dose_allocation;                   // what proportion of total daily doses are reserved for reactive strategies

        std::vector<int> min_age;
        // double revaccination_follow_up;                              // what proportion of people follow-up for second doses

        // void _geographic_cluster_vaccination();
};


#endif
