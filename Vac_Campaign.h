// Vac_Campaign.h
#ifndef __VAC_CAMPAIGN_H
#define __VAC_CAMPAIGN_H

#include "Person.h"
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
        virtual ~Vac_Campaign() {} //is this correct to address a vtable error?

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

        size_t get_urgent_queue_size() { return urgent_queue.size(); }
        size_t get_standard_queue_size() { return standard_queue.size(); }
        size_t get_revaccinate_queue_size(int day) { return revaccinate_queue[day].size(); }

        void prioritize_vaccination(Person* p) { urgent_queue.push(p); }
        void schedule_vaccination(Person* p) { standard_queue.push(p); }

        int get_doses_available(int day, VaccineAllocationType alloc) { return doses_available[day][alloc]; }
        void set_doses_available( std::vector< std::vector<int> > da ) {
            doses_available = da;
            //doses_used.resize(da.size(), std::vector<int>(NUM_OF_VACCINATION_QUEUE_TYPES));
            revaccinate_queue.resize(da.size());
            queue_tally.resize(da.size(), vector<int>(NUM_OF_VACCINATION_QUEUE_TYPES));
        }

        std::vector<int> get_min_age() const  { return min_age; }
        int get_min_age(size_t day)    const  { return min_age.at(day); }
        void set_min_age(std::vector<int> ma) { min_age = ma; }

        bool is_age_eligible_on(int age, int day) { return age >= min_age.at(day); }

        Vaccinee* next_vaccinee(size_t day) {
            // Set up defaults
            Person* person = nullptr;
            VaccineAllocationType source_of_dose = NUM_OF_VACCINE_ALLOCATION_TYPES;
            VaccinationQueueType  vaccinee_queue = NUM_OF_VACCINATION_QUEUE_TYPES;

            // Do we have doses available today?
            const size_t urgent_doses_avail      = unlim_urgent_doses ? 1 : doses_available.at(day)[URGENT_ALLOCATION];
            const size_t standard_doses_avail    = doses_available.at(day)[STANDARD_ALLOCATION];

            if (urgent_doses_avail + standard_doses_avail == 0) { return nullptr; } // might as well bail

            // How many in each queue are waiting to be vaccinated today?
            const size_t urgent_size = urgent_queue.size();
            const size_t std_size = standard_queue.size();
            const size_t revac_size = revaccinate_queue[day].size();

            // TODO -- we may want to have different priorities for finishing a multi-dose vaccine course versus
            // administering (e.g. annual) booster doses
            // Possible TODO -- this does not currently handle the scenario where the priority is urgent, then
            // revaccinations, then standard vaccinations.
            if (revac_size and (not prioritize_first_doses or not (urgent_size or std_size))) {
            // if (revac_size and (not prioritize_first_doses and not (urgent_size or std_size))) {
                // we have someone to revaccinate, and either we're prioritizing completing courses
                // or we have no one neededing a first dose
                source_of_dose = standard_doses_avail ? STANDARD_ALLOCATION
                                 : flexible_queue_allocation and urgent_doses_avail ? URGENT_ALLOCATION
                                 : NUM_OF_VACCINE_ALLOCATION_TYPES;

                if (source_of_dose < NUM_OF_VACCINE_ALLOCATION_TYPES) { // i.e., urgent or standard
                    person = *(revaccinate_queue[day].begin());
                    revaccinate_queue[day].erase(person);
                    vaccinee_queue = REVACCINATE_QUEUE;
                }
            }

            if (urgent_size and not person) {
                source_of_dose = urgent_doses_avail ? URGENT_ALLOCATION
                                 : flexible_queue_allocation and standard_doses_avail ? STANDARD_ALLOCATION
                                 : NUM_OF_VACCINE_ALLOCATION_TYPES;

                if (source_of_dose < NUM_OF_VACCINE_ALLOCATION_TYPES) {
                    person = urgent_queue.front();
                    urgent_queue.pop();
                    vaccinee_queue = URGENT_QUEUE;
                }
                if(unlim_urgent_doses) { doses_available.at(day)[URGENT_ALLOCATION]++; }
            }

            if (std_size and not person) {
                source_of_dose = standard_doses_avail ? STANDARD_ALLOCATION
                                 : flexible_queue_allocation and urgent_doses_avail ? URGENT_ALLOCATION
                                 : NUM_OF_VACCINE_ALLOCATION_TYPES;

                if (source_of_dose < NUM_OF_VACCINE_ALLOCATION_TYPES) {
                    person = standard_queue.front();
                    standard_queue.pop();
                    vaccinee_queue = STANDARD_QUEUE;
                }
            }

            // these lines could be combined, but I think the intent is clearer this way
            Vaccinee* vaccinee = person ? new Vaccinee(person, vaccinee_queue, source_of_dose) : nullptr;
            return vaccinee;
        }

        bool vaccinate(Vaccinee* v, int day) { return is_age_eligible_on(v->get_person()->getAge(), day) and v->vaccinate(day); }

        void tally_dose(size_t day, Vaccinee* vac) {
            doses_available.at(day)[vac->dose_source]--;
            queue_tally[day][vac->status]++;
        }
        int get_dose_tally(size_t day, VaccinationQueueType q) { return queue_tally[day][q]; }

        void schedule_revaccination(size_t day, Vaccinee* vac) {
            if (revaccinate_queue.size() > day) {
                revaccinate_queue[day].insert(vac->person);
                if(vac->get_status() == URGENT_QUEUE and unlim_urgent_doses) {
                    doses_available.at(day)[STANDARD_ALLOCATION]++;
                }
            }
        }

        void reschedule_remaining_revaccinations(size_t today) {
            const size_t queue_length_in_days = revaccinate_queue.size(); // how many days long is the simulation
            std::set<Person*> remaining_people = revaccinate_queue[today];
            if (remaining_people.size() and queue_length_in_days > today + 1) {
                revaccinate_queue[today + 1].insert(remaining_people.begin(), remaining_people.end());
                revaccinate_queue[today].clear();
            }
        }

        void ring_scheduling(int day, vector< set<Person*> > tracedContacts);
        void geographic_scheduling(int day, vector< set<Person*> > targetedPeople, Community* community);
        void location_scheduling(int day, vector< set<Person*> > targetedPeople);

        void reactive_strategy(int day, vector< set<Person*> > targetedPeople, Community* community);

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
        std::vector<int> start_of_campaign;                    // index by VacCampaignType and returns start date as sim day
        std::vector<int> end_of_campaign;                      // index by VacCampaignType and returns start date as sim day

        std::queue<Person*> urgent_queue;                      // people who need to be vaccinated, determined during simulation
        std::queue<Person*> standard_queue;                    // people who will be vaccinated, known before transmission sim begins
        std::vector< std::set<Person*> > revaccinate_queue;    // indexed by sim day

        bool prioritize_first_doses;                           // do we prioritize first doses, or completing vac schedules?
        bool flexible_queue_allocation;                        // can doses be used from any allocation, or only as intended?
        bool unlim_urgent_doses;                               // should an unlimited number of doses be used for the urgent queue (i.e. reactive strategy)?

        std::vector< std::vector<int> > doses_available;       // for each day, number of doses available for each queue
        std::vector< std::vector<int> > queue_tally;           // for each day, number of doses used for each queue

        VacCampaignType reactive_vac_strategy;                 // parameter for type of reactive strategy (if one is active)
        double reactive_vac_dose_allocation;                   // what proportion of total daily doses are reserved for reactive strategies

        std::vector<int> min_age;
        // double revaccination_follow_up;                              // what proportion of people follow-up for second doses

        // void _geographic_cluster_vaccination();
};


#endif
