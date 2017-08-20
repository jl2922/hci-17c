#include "det.h"

bool operator==(const Det& lhs, const Det& rhs) { return lhs.up == rhs.up && lhs.dn == rhs.dn; }