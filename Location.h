// Location.h
// A location, which is basically a building.
// Maintains lists of which humans are in the building at 3 points in the day
#ifndef __LOCATION_H
#define __LOCATION_H

#include <queue>

class Person;

class Location {
    public:
        Location();
        virtual ~Location();
        int getID() const { return _ID; }
        void setType(LocationType t) { _type = t; }
        void setEssential(bool e) { _essential = e; }
        bool isEssential() const { return _essential; }
        LocationType getType() const { return _type; }
        void addPerson(Person *p) { _person.push_back(p); }
        bool removePerson(Person *p);
        int getNumPeople() const { return _person.size(); }
        std::vector<Person*> getPeople() { return _person; }
//        Person* findMom();                                            // Try to find a resident female of reproductive age
        void addNeighbor(Location *p);
        int getNumNeighbors() const { return _neighbors.size(); }
//        Location *getNeighbor(int n) { return _neighbors[n]; }
        inline Person* getPerson(int idx) { return _person[idx]; }
        void setCoordinates(std::pair<double, double> c) { _coord = c; }
        std::pair<double, double> getCoordinates() { return _coord; }
        void setX(double x) { _coord.first = x; }
        void setY(double y) { _coord.second = y; }
        double getX() const { return _coord.first; }
        double getY() const { return _coord.second; }

        bool operator == ( const Location* other ) const { return _ID == other->_ID; }

        // We use this to make sure that locations are iterated through in a well-defined order (by ID), rather than by mem address
        struct LocPtrComp { bool operator()(const Location* A, const Location* B) const { return A->getID() < B->getID(); } };

    protected:
        int _ID;                                                      // original identifier in location file
        bool _essential;
        LocationType _type;
        std::vector<Person*> _person;                                 // pointers to people who come to this location
        std::set<Location*, Location::LocPtrComp> _neighbors;
        static size_t NEXT_ID;                                        // unique ID to assign to the next Location allocated
        std::pair<double, double> _coord;                             // (x,y) coordinates for location
};



#endif
