#ifndef HCI_WAVEFUNCTION_H_
#define HCI_WAVEFUNCTION_H_

#include "../std.h"

#include "det.h"
#include "term.h"

class Wavefunction {
 private:
  std::list<Term> terms;

 public:
  Wavefunction() {}

  size_t size() { return terms.size(); }

  void append_term(const Det& det, const double coef) { terms.push_back(Term(det, coef)); }

  const std::list<Term>& get_terms() const { return terms; }

  void set_coefs(const std::vector<double>& coefs) {
    size_t i = 0;
    for (auto& term : terms) term.coef = coefs[i++];
  }

  std::vector<Det> get_dets() const {
    std::vector<Det> dets;
    dets.reserve(terms.size());
    for (const auto& term : terms) dets.push_back(term.det);
    return dets;
  }

  std::vector<double> get_coefs() const {
    std::vector<double> coefs;
    coefs.reserve(terms.size());
    for (const auto& term : terms) coefs.push_back(term.coef);
    return coefs;
  }

  void sort_by_coefs() {
    terms.sort([](const Term& a, const Term& b) -> bool { return fabs(a.coef) > fabs(b.coef); });
  }

  void clear() { terms.clear(); }
};

#endif