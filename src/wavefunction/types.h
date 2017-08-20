#ifndef TYPES_H_
#define TYPES_H_

#include "../std.h"

typedef std::uint16_t Orbital;  // [0, 65535 (2^16-1)].
typedef std::vector<Orbital> Orbitals;
typedef std::pair<Orbital, Orbital> OrbitalPair;
typedef std::pair<Orbitals, Orbitals> OrbitalsPair;

#endif