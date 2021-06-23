// Vac_Campaign.h
#ifndef __VAC_CAMPAIGN_H
#define __VAC_CAMPAIGN_H

#include "Person.h"
#include <queue>

class Person;

enum VaccinationQueueType {
    URGENT_QUEUE;
    STANDARD_QUEUE;
    REVACCINATE_QUEUE;
    NUM_OF_VACCINATION_QUEUE_TYPES;
};

enum VaccineAllocationType {
    URGENT_ALLOCATION;
    STANDARD_ALLOCATION;
    NUM_OF_VACCINE_ALLOCATION_TYPES;
};

class Vaccinee {
    public:
        Vaccinee() : person(nullptr), status(NUM_OF_VACCINATION_QUEUE_TYPES), dose_source(NUM_OF_VACCINE_ALLOCATION_TYPES) {}
        Vaccinee(Person* p, VaccinationQueueType stat, VaccineAllocationType source) : person(p), status(stat), dose_source(source) {}
        bool vaccinate (int day) { return person->vaccinate(day); } // true if person was vaccinatable
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
};

class Vac_Campaign {
    public:
        Vac_Campaign() {
            prioritize_first_doses = false;
            flexible_queue_allocation = false;
        }
        virtual ~Vac_Campaign();

        void set_prioritize_first_doses(bool val) { prioritize_first_doses = val; }
        void set_flexible_queue_allocation(bool val) { flexible_queue_allocation = val; }

        void prioritize_vaccination(Person* p) { urgent_queue.push(p); }
        void schedule_vaccination(Person* p) { standard_queue.push(p); }

        void set_doses_available( std::vector< std::vector<int> > da ) {
            doses_available = da;
            doses_used.resize(da.size(), std::vector<int>(NUM_OF_VACCINATION_QUEUE_TYPES));
            revaccinate_queue.resize(da.size());
        }

        Vaccinee* next_vaccinee(size_t day) {
            // Set up defaults
            Person* person = nullptr;
            VaccineAllocationType source_of_dose = NUM_OF_VACCINE_ALLOCATION_TYPES;
            VaccinationQueueType  vaccinee_queue = NUM_OF_VACCINATION_QUEUE_TYPES;

            // Do we have doses available today?
            const size_t urgent_doses_avail      = doses_available.at(day)[URGENT_ALLOCATION];
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
                // we have someone to revaccinate, and either we're prioritizing completing courses
                // or we have no one neededing a first dose
                source_of_dose = standard_doses_avail ? STANDARD_ALLOCATION
                                 : flexible_queue_allocation and urgent_doses_avail ? URGENT_ALLOCATION
                                 : NUM_OF_VACCINE_ALLOCATION_TYPES;

                if (source_of_dose < NUM_OF_VACCINE_ALLOCATION_TYPES) { // i.e., urgent or standard
                    person = *(revaccinate_queue[day].begin());
                    revaccinate_queue[day].erase(p);
                    vaccinee_queue = REVACCINATE_QUEUE;
                }
            } else if (urgent_queue.size()) {
                source_of_dose = urgent_doses_avail ? URGENT_ALLOCATION
                                 : flexible_queue_allocation and standard_doses_avail ? STANDARD_ALLOCATION
                                 : NUM_OF_VACCINE_ALLOCATION_TYPES;

                if (source_of_dose < NUM_OF_VACCINE_ALLOCATION_TYPES) {
                    person = urgent_queue.front();
                    urgent_queue.pop();
                    vaccinee_queue = URGENT_QUEUE;
                }
            } else if (standard_queue.size()) {
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

        void tally_dose(size_t day, Vaccinee* vac) {
            doses_available.at(day)[vac->dose_source]--;
            queue_tally[day][vac->status]++;
        }

        void schedule_revaccination(size_t day, Vaccinee* vac) {
            revaccinate_queue.at(day).add(vac->person);
        }

        void reschedule_remaining_revaccinations(size_t day) {
            const size_t queue_length_in_days = revaccinate_queue.size();
            std::set<Person*> remaining_people = revaccinate_queue[day];
            if (remaining_people.size() and queue_length_in_days > day + 1) {
                revaccinate_queue[day + 1].insert(remaining_queue);
                revaccinate_queue[day].clear();
            }
        }

    private:
        std::queue<Person*> urgent_queue;                      // people who need to be vaccinated, determined during simulation
        std::queue<Person*> standard_queue;                    // people who will be vaccinated, known before transmission sim begins
        std::vector< std::set<Person*> > revaccinate_queue;    // indexed by sim day
        bool prioritize_first_doses;                           // do we prioritize first doses, or completing vac schedules?
        std::vector< std::vector<int> > doses_available;       // for each day, number of doses available for each queue
        std::vector< std::vector<int> > queue_tally;           // for each day, number of doses used for each queue
        bool flexible_queue_allocation;                        // can doses be used from any allocation, or only as intended?
};


#endif
