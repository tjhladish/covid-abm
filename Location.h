// Location.h
// A location, which is basically a building.
// Maintains lists of which humans are in the building at 3 points in the day
#ifndef __LOCATION_H
#define __LOCATION_H

#include <queue>

class Person;

class Location {
    // TODO -- create derived classes for different types of locations
    public:
        Location();
        Location(const Location& o) {
            _ID                       = o._ID;
            _essential                = o._essential;
            _riskiness                = o._riskiness;
            _public_transmission_risk = o._public_transmission_risk;
            _type                     = o._type;
            _person                   = o._person;
            _visitors                 = o._visitors;
            _visit_durations          = o._visit_durations;
            _neighbors                = o._neighbors;
            _hospital                 = o._hospital;
            _coord                    = o._coord;
        };

        virtual ~Location();

        struct LocPtrComp { bool operator()(const Location* A, const Location* B) const { return A->getID() < B->getID(); } };

        int getID() const { return _ID; }
        static void reset_ID_counter() { NEXT_ID = 0; } // for use by community destructor

        void setType(LocationType t) { _type = t; }
        void setEssential(bool e) { _essential = e; }
        bool isEssential() const { return _essential; }
        bool isNonEssential() const { return !_essential; }
        LocationType getType() const { return _type; }
        void addPerson(Person *p) { _person.push_back(p); }
        void addVisitor(Person *p, double dur) { _visitors.push_back(p); _visit_durations.push_back(dur); }
        void clearVisitors() { _visitors.clear(); _visit_durations.clear(); }                     // likey don't want to resize
        bool removePerson(Person *p);
        int getNumPeople() const { return _person.size(); }             // employees, residents, etc.; not including visitors
        int getNumVisitors() const { return _visitors.size(); }          // visitors change daily
        std::vector<Person*> getPeople() { return _person; }
        void setPeople(std::vector<Person*> p) { _person = p; }
        std::vector<Person*> getVisitors() { return _visitors; }
        void setVisitors(std::vector<Person*> v) { _visitors = v; }
        std::vector<double> getVisitDurations() { return _visit_durations; }
//        Person* findMom();                                            // Try to find a resident female of reproductive age
        void addNeighbor(Location* loc);
        int getNumNeighbors() const { return _neighbors.size(); }
        std::set<Location*, LocPtrComp> getNeighbors() { return _neighbors; }
        void setNeighbors(std::set<Location*, LocPtrComp> n) { _neighbors = n; }
        void setHospital(Location* hosp) { _hospital = hosp; }
        Location* getHospital() const { return _hospital; }
        inline Person* getPerson(int idx) { return _person[idx]; }
        void setCoordinates(std::pair<double, double> c) { _coord = c; }
        std::pair<double, double> getCoordinates() { return _coord; }
        void setX(double x) { _coord.first = x; }
        void setY(double y) { _coord.second = y; }
        double getX() const { return _coord.first; }
        double getY() const { return _coord.second; }
        void setRiskiness(float ra) { _riskiness = ra; }
        float getRiskiness() const { return _riskiness; }
        void setPublicTransmissionRisk(PublicTransmissionType pt) { _public_transmission_risk = pt; }
        PublicTransmissionType getPublicTransmissionRisk() const { return _public_transmission_risk; }

        bool operator == ( const Location* other ) const { return _ID == other->_ID; }
        void dumper() const;

        // We use this to make sure that locations are iterated through in a well-defined order (by ID), rather than by mem address

    protected:
        // TODO -- figure out how to represent people coming to a place for different reasons,
        // e.g. hospital patients/staff, nursinghome residents/staff
        int _ID;                                                      // original identifier in location file
        bool _essential;
        float _riskiness;                                             // score on U(0,1) of household's risk threshold
        PublicTransmissionType _public_transmission_risk;             // is this a location where there is (high, low, or no) risk of tranmission to e.g. customers
        LocationType _type;
        std::vector<Person*> _person;                                 // people who predictably come to this location
        std::vector<Person*> _visitors;                               // people who probabilistically visit this location (currently only for businesses)
        std::vector<double> _visit_durations;                         // people who probabilistically visit this location (currently only for businesses)
        std::set<Location*, LocPtrComp> _neighbors;
        Location* _hospital;                                          // for houses, the associated hospital (nullptr for others)
        static size_t NEXT_ID;                                        // unique ID to assign to the next Location allocated
        std::pair<double, double> _coord;                             // (x,y) coordinates for location
};



#endif
