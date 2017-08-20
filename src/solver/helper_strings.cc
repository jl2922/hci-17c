#include "helper_strings.h"

#include "../config.h"
#include "../time.h"
#include "omp.h"

void HelperStrings::setup(const std::vector<Det>& dets) {
  this->dets = dets;
  setup_ab();
  setup_ab_m1();

  size_t n_dets = dets.size();
  cached.assign(n_dets, false);
  cached_connections.resize(n_dets);
  cache_size = Config::get<size_t>("cache_size", 1000);

#pragma omp parallel
  {
    const int thread_id = omp_get_thread_num();
    if (thread_id == 0) {
      connected.resize(omp_get_num_threads());
      one_up.resize(omp_get_num_threads());
    }
#pragma omp barrier
    connected[thread_id].assign(n_dets, false);
    one_up[thread_id].assign(n_dets, false);
  }
}

void HelperStrings::setup_ab() {
  for (std::size_t i = 0; i < dets.size(); i++) {
    ab[dets[i].up.encode()].first.push_back(i);
    ab[dets[i].dn.encode()].second.push_back(i);
  }
}

void HelperStrings::setup_ab_m1() {
  for (std::size_t i = 0; i < dets.size(); i++) {
    const auto& up_elecs = dets[i].up.get_elec_orbs();
    SpinDet det_up(dets[i].up);
    for (std::size_t j = 0; j < up_elecs.size(); j++) {
      det_up.set_orb(up_elecs[j], false);
      ab_m1[det_up.encode()].first.push_back(i);
      det_up.set_orb(up_elecs[j], true);
    }

    const auto& dn_elecs = dets[i].dn.get_elec_orbs();
    SpinDet det_dn(dets[i].dn);
    for (std::size_t j = 0; j < dn_elecs.size(); j++) {
      det_dn.set_orb(dn_elecs[j], false);
      ab_m1[det_dn.encode()].second.push_back(i);
      det_dn.set_orb(dn_elecs[j], true);
    }
  }
}

std::vector<std::pair<size_t, double>> HelperStrings::find_connections(const std::size_t i) {
  if (cached[i]) return cached_connections[i];

  const int thread_id = omp_get_thread_num();
  std::vector<std::pair<size_t, double>> connections;
  const Det& det = dets[i];
  const auto& up_elecs = det.up.get_elec_orbs();
  const auto& dn_elecs = det.dn.get_elec_orbs();

  // Two up/dn excitations.
  if (ab.find(det.dn.encode()) != ab.end()) {
    for (const std::size_t det_id : ab.find(det.dn.encode())->second.second) {
      if (!connected[thread_id][det_id]) {
        if (det_id < i) continue;
        connected[thread_id][det_id] = true;
        const double H = hamiltonian(det, dets[det_id]);
        connections.push_back(std::make_pair(det_id, H));
      }
    }
  }
  if (ab.find(det.up.encode()) != ab.end()) {
    for (const std::size_t det_id : ab.find(det.up.encode())->second.first) {
      if (!connected[thread_id][det_id]) {
        if (det_id < i) continue;
        connected[thread_id][det_id] = true;
        const double H = hamiltonian(det, dets[det_id]);
        connections.push_back(std::make_pair(det_id, H));
      }
    }
  }

  // One up one dn excitation.
  SpinDet det_up(det.up);
  SpinDet det_dn(det.dn);
  std::vector<std::size_t> one_ups;

  for (std::size_t k = 0; k < up_elecs.size(); k++) {
    det_up.set_orb(up_elecs[k], false);
    const auto& key_up = det_up.encode();
    const auto& kv_up = ab_m1.find(key_up);
    if (kv_up != ab_m1.end()) {
      for (const std::size_t det_id : kv_up->second.first) {
        if (det_id < i) continue;
        one_up[thread_id][det_id] = true;
        one_ups.push_back(det_id);
      }
    }
    det_up.set_orb(up_elecs[k], true);
  }

  for (std::size_t k = 0; k < dn_elecs.size(); k++) {
    det_dn.set_orb(dn_elecs[k], false);
    const auto& key_dn = det_dn.encode();
    if (ab_m1.find(key_dn) != ab_m1.end()) {
      for (const std::size_t det_id : ab_m1.find(key_dn)->second.second) {
        if (one_up[thread_id][det_id] && !connected[thread_id][det_id]) {
          connected[thread_id][det_id] = true;
          const double H = hamiltonian(det, dets[det_id]);
          connections.push_back(std::make_pair(det_id, H));
        }
      }
    }
    det_dn.set_orb(dn_elecs[k], true);
  }

  // Reset connected and return.
  for (const std::size_t det_id : one_ups) one_up[thread_id][det_id] = false;
  for (const auto& connection : connections) connected[thread_id][connection.first] = false;

  if (connections.size() < cache_size) {
    cached_connections[i] = connections;
    cached[i] = true;
  }

  return connections;
}
