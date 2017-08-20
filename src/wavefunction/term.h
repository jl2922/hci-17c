#ifndef TERM_H_
#define TERM_H_

#include "../std.h"
#include "det.h"

class Term {
 public:
  Det det;
  double coef;

  Term(const Det& det, const double coef) {
    this->det = det;
    this->coef = coef;
  }
};

#endif