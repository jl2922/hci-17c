#include "solver.h"

#include <boost/functional/hash.hpp>
#include <boost/format.hpp>
#include "../parallel.h"
#include "../std.h"
#include "../time.h"
#include "../wavefunction/wavefunction.h"
#include "davidson.h"

Det Solver::generate_hf_det() {
  Det det;
  for (size_t i = 0; i < n_up; i++) det.up.set_orb(i, true);
  for (size_t i = 0; i < n_dn; i++) det.dn.set_orb(i, true);
  return det;
}

void Solver::variation(const double eps_var) {
  const double THRESHOLD = 1.0e-6;

  // Setup HF or existing wf as initial wf and evaluate energy.
  if (wf.size() == 0) {
    const Det& det_hf = generate_hf_det();
    wf.append_term(det_hf, 1.0);
    energy_hf = energy_var = hamiltonian(det_hf, det_hf);
    if (Parallel::is_master()) printf("HF energy: %#.15g Ha\n", energy_hf);
  }

  std::unordered_set<OrbitalsPair, boost::hash<OrbitalsPair>> var_dets_set;
  std::list<Det> new_dets;
  for (const auto& term : wf.get_terms()) var_dets_set.insert(term.det.encode());
  double energy_var_new = 0.0;  // Ensures the first iteration will run.
  size_t n_iter = 0;
  converged = false;
  while (!converged) {
    Time::start("Variation: " + std::to_string(n_iter));

    for (const auto& term : wf.get_terms()) {
      const auto& connected_dets = find_connected_dets(term.det, eps_var / fabs(term.coef));
      for (const auto& new_det : connected_dets) {
        const auto new_det_code = new_det.encode();
        if (var_dets_set.count(new_det_code) == 0) {
          var_dets_set.insert(new_det_code);
          new_dets.push_back(new_det);
        }
      }
    }

    if (Parallel::is_master()) {
      printf("New / total dets: %'zu / %'zu\n", new_dets.size(), var_dets_set.size());
    }
    Time::checkpoint("found new dets");

    for (const auto& new_det : new_dets) {
      wf.append_term(new_det, 0.0);
    }

    energy_var_new = diagonalize(new_dets.size() > 0);
    if (fabs(energy_var - energy_var_new) < THRESHOLD) converged = true;
    energy_var = energy_var_new;
    if (Parallel::is_master()) {
      printf("Variation energy: %#.15g Ha\n", energy_var);
      printf("Correlation energy: %#.15g Ha\n", energy_var - energy_hf);
    }

    new_dets.clear();

    Time::end();
    n_iter++;
  }

  if (Parallel::is_master()) {
    printf("Final variation energy: %#.15g Ha\n", energy_var);
    printf("Correlation energy: %#.15g Ha\n", energy_var - energy_hf);
  }
}

double Solver::diagonalize(const bool has_new_dets) {
  std::vector<double> diagonal;
  std::vector<double> initial_vector;
  const size_t max_iterations = has_new_dets ? 5 : 10;
  diagonal.reserve(wf.size());
  initial_vector.reserve(wf.size());
  for (const auto& term : wf.get_terms()) {
    const auto& det = term.det;
    diagonal.push_back(hamiltonian(det, det));
    initial_vector.push_back(term.coef);
  }

  std::function<double(Det, Det)> hamiltonian_func =
      std::bind(&Solver::hamiltonian, this, std::placeholders::_1, std::placeholders::_2);
  HelperStrings helper_strings(hamiltonian_func);
  helper_strings.setup(wf.get_dets());
  Time::checkpoint("helper strings generated");
  std::function<std::vector<double>(std::vector<double>)> apply_hamiltonian_func =
      std::bind(&Solver::apply_hamiltonian, this, std::placeholders::_1, helper_strings);

  Davidson davidson(diagonal, apply_hamiltonian_func, wf.size());
  if (Parallel::is_master()) davidson.set_verbose(true);
  const size_t n_iter = davidson.diagonalize(initial_vector, max_iterations);
  if (!has_new_dets && n_iter < max_iterations) converged = true;

  const double energy_var = davidson.get_lowest_eigenvalue();
  const auto& coefs_new = davidson.get_lowest_eigenvector();

  wf.set_coefs(coefs_new);
  wf.sort_by_coefs();

  return energy_var;
}

#pragma omp declare reduction(      \
    vec_double_plus : std::vector < \
    double > : std::transform(      \
                 omp_out            \
                     .begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus < double > ())) initializer(omp_priv = omp_orig)

std::vector<double> Solver::apply_hamiltonian(
    const std::vector<double>& vec, HelperStrings& helper_strings) {
  const std::size_t n_dets = vec.size();
  const size_t proc_id = Parallel::get_id();
  const size_t n_procs = Parallel::get_n();
  std::vector<double> res(n_dets, 0.0);

  const auto& dets = wf.get_dets();

#pragma omp parallel for reduction(vec_double_plus : res) schedule(dynamic, 10)
  for (size_t i = proc_id; i < n_dets; i += n_procs) {
    const auto& connections = helper_strings.find_connections(i);
    for (const auto connection : connections) {
      const size_t j = connection.first;
      const double H_ij = connection.second;
      res[i] += H_ij * vec[j];
      if (i != j) {
        res[j] += H_ij * vec[i];
      }
    }
  }

  Parallel::reduce_to_sum_vector(res);
  Time::checkpoint("hamiltonian applied");

  return res;
}


void Solver::save_variation_result(const std::string& filename) {
  if (Parallel::is_master()) {
    std::ofstream var_file;
    var_file.open(filename);
    var_file << boost::format("%.17g %.17g\n") % energy_hf % energy_var;
    var_file << boost::format("%d %d %d\n") % n_up % n_dn % wf.size();
    for (const auto& term : wf.get_terms()) {
      var_file << boost::format("%.17g\n") % term.coef;
      var_file << term.det.up << std::endl << term.det.dn << std::endl;
    }
    var_file.close();
    printf("Variation result saved to: %s\n", filename.c_str());
  }
}


bool Solver::load_variation_result(const std::string& filename) {
  std::ifstream var_file;
  size_t n_dets;
  Orbital orb_id;
  double coef;
  var_file.open(filename);
  if (!var_file.is_open()) return false;  // Does not exist.
  var_file >> energy_hf >> energy_var;
  var_file >> n_up >> n_dn >> n_dets;
  wf.clear();
  for (std::size_t i = 0; i < n_dets; i++) {
    var_file >> coef;
    Det det;
    for (std::size_t j = 0; j < n_up; j++) {
      var_file >> orb_id;
      det.up.set_orb(orb_id, true);
    }
    for (std::size_t j = 0; j < n_dn; j++) {
      var_file >> orb_id;
      det.dn.set_orb(orb_id, true);
    }
    wf.append_term(det, coef);
  }
  var_file.close();
  if (Parallel::is_master()) {
    printf("Loaded %'zu dets from: %s\n", n_dets, filename.c_str());
    printf("HF energy: %#.15g Ha\n", energy_hf);
    printf("Variation energy: %#.15g Ha\n", energy_var);
    printf("Correlation energy: %#.15g Ha\n", energy_var - energy_hf);
  }
  return true;
}