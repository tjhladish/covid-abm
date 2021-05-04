// Location.cpp

#include <cstdlib>
#include <cstring>
#include <climits>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Person.h"
#include "Location.h"
#include "Parameters.h"

using namespace covid::standard;

size_t Location::NEXT_ID = 0;

Location::Location()
    : _person(0) {
    _ID = NEXT_ID++;
    _coord = make_pair(0.0, 0.0);
    _type = NUM_OF_LOCATION_TYPES; // compileable, but not sensible value, because it must be set elsewhere
    _essential = true;
    _riskiness= 0.0;
    _hospital = nullptr;
}


Location::~Location() {
    _person.clear();
    _neighbors.clear();
}


void Location::dumper() const {
    cerr << "Location ID: " << _ID << endl;
    cerr << "\ttype: " << _type << endl;
    cerr << "\tcoords: " << _coord.first << ", " << _coord.second << endl;
    cerr << "\tessential: " << _essential << endl;
    cerr << "\triskiness: " << _riskiness << endl;
    cerr << "\tpeople: "; /*for (auto p: _person) { cerr << p->getID() << " "; }*/ cerr << _person.size(); cerr << endl;
    cerr << "\tvisitors: "; /*for (auto p: _visitors) { cerr << p->getID() << " "; }*/ cerr << _visitors.size(); cerr << endl;
    cerr << "\tneighbors: "; for (auto l: _neighbors) { cerr << l->getID() << " "; } cerr << endl;
}


bool Location::removePerson(Person* p) { // TODO -- consider implications of using set for container, to make removing faster
    for (unsigned int i=0; i<_person.size(); i++) {
        if (_person[i] == p) {
            _person[i] = _person.back();
            _person.pop_back();
            return true;
        }
    }
    return false;
}

/*
// Calling scope must verify that returned person is not nullptr
Person* Location::findMom() {
    vector<Person*> residents = getResidents();
    vector<Person*> potential_moms;
    int minage = 15;
    int maxage = 45;
    for (auto p: residents) {
        if (p->getSex() == FEMALE and p->getAge() >= minage and p->getAge() <= maxage) potential_moms.push_back(p);
    }
    if (potential_moms.size() == 0) return nullptr;
    int r = gsl_rng_uniform_int(RNG, potential_moms.size());
    Person* mom = potential_moms[r];
    return mom;
}*/


// addNeighbor - adds location p to the location's neighbor list.
// Note that this relationship is one-way.
void Location::addNeighbor(Location* loc) {
    _neighbors.insert(loc);
}
