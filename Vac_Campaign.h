// Vac_Campaign.h
#ifndef __VAC_CAMPAIGN_H
#define __VAC_CAMPAIGN_H

#include <queue>

class Person;

enum VaccinationQueueType {
    URGENT_QUEUE;
    STANDARD_QUEUE;
    NUM_OF_VACCINATION_QUEUE_TYPES;
};

enum VaccineAllocationType {
    URGENT_ALLOCATION;
    STANDARD_ALLOCATION;
    NUM_OF_VACCINE_ALLOCATION_TYPES;
};

struct Vaccinee {
    Vaccinee() : person(nullptr), status(NUM_OF_VACCINATION_QUEUE_TYPES), dose_source(NUM_OF_VACCINE_ALLOCATION_TYPES) {}
    Vaccinee(std::queue<Person*> q, VaccinationQueueType stat, VaccineAllocationType source) : person(q.front()), status(stat), dose_source(source) { q.pop(); }
    Person* person;
    // We allow the possibility that the dose source won't match the vaccinee status because we may want to
    // consider idealized scenarios where doses can be rapidly reallocated between vaccine distribution
    // programs.
    VaccinationQueueType status;       // was this vaccinee standard (planned, proactive) or urgent (reactive, traced, etc.)
    VaccineAllocationType dose_source; // was this dose allocated to the urgent or standard queue
};

class Vac_Campaign {
    public:
        Vac_Campaign();
        virtual ~Vac_Campaign();

        void prioritize_vaccination(Person* p) { urgent_queue.push(p); }
        void schedule_vaccination(Person* p) { standard_queue.push(p); }

        Vaccinee* next_vaccinee(size_t day) {
            Vaccinee* vaccinee = nullptr;//new Vaccinee();
CALL DELETE VACCINEE WHEN WE DETERMINE WHETHER TO VACCINATE THIS PERSON
            VaccineAllocationType source_of_dose = NUM_OF_VACCINE_ALLOCATION_TYPES;
            const size_t urgent_doses_avail  = doses_available[day][URGENT_ALLOCATION];
            const size_t standard_doses_avail = doses_available[day][STANDARD_ALLOCATION];
            if (urgent_queue.size()) {
                if (urgent_doses_avail) {
                    vaccinee = new Vaccinee(urgent_queue, URGENT_QUEUE, URGENT_ALLOCATION);
                } else if (not strict_queue_allocation and standard_doses_avail) {
                    vaccinee = new Vaccinee(urgent_queue, URGENT_QUEUE, STANDARD_ALLOCATION);
                }
            } else if (standard_queue.size()) {
                if (standard_doses_avail) {
                    vaccinee = new Vaccinee(standard_queue, STANDARD_QUEUE, STANDARD_ALLOCATION);
                } else if (not strict_queue_allocation and urgent_doses_avail) {
                    vaccinee = new Vaccinee(standard_queue, STANDARD_QUEUE, URGENT_ALLOCATION);
                }
            }
            return vaccinee;
        }

        void tally_dose(size_t day, Vaccinee* vac) {
            doses_available[day][vac.dose_source]--;
            doses_used[day][vac.status]++;
        }
        

    private:
        std::queue<Person*> urgent_queue;                   // people who need to be vaccinated, determined during simulation
        std::queue<Person*> standard_queue;                 // people who will be vaccinated, known before transmission sim begins
        bool priortize_fist_doses;                          // do we prioritize first doses, or completing vac schedules?
        std::vector< vector<int> > doses_available;         // for each day, number of doses available for each queue
        std::vector< vector<int> > doses_used;              // for each day, number of doses used for each queue
        bool strict_queue_allocation;                       // can doses only be used for their allocated queue?
};


#endif
