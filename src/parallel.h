#ifndef PARALLEL_H_
#define PARALLEL_H_

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "omp.h"
#include "std.h"

#ifndef SERIAL

template <class T>
class vector_plus : public std::binary_function<std::vector<T>, std::vector<T>, std::vector<T>> {
 public:
  std::vector<T> operator()(const std::vector<T>& lhs, const std::vector<T>& rhs) const {
    std::vector<T> v(lhs.size());
    std::transform(lhs.begin(), lhs.end(), rhs.begin(), v.begin(), std::plus<T>());
    return (v);
  }
};

class Parallel {
 private:
  size_t id;
  size_t n;
  boost::mpi::environment* env;  // For MPI 1.1.
  boost::mpi::communicator world;

  Parallel() {
    id = world.rank();
    n = world.size();
  }

  // Singleton pattern boilerplate.
  static Parallel& get_instance() {
    static Parallel instance;
    return instance;
  }

 public:
  static void init(boost::mpi::environment& env) { Parallel::get_instance().env = &env; }

  static size_t get_id() { return Parallel::get_instance().id; }

  static size_t get_n() { return Parallel::get_instance().n; }

  static bool is_master() { return Parallel::get_instance().id == 0; }

  static std::string get_host() { return Parallel::get_instance().env->processor_name(); }

  static void barrier() {
    fflush(stdout);
    Parallel::get_instance().world.barrier();
  }

  static void print_info() {
    Parallel::barrier();
#pragma omp parallel
    {
      if (omp_get_thread_num() == 0) {
        printf(
            "Proc %zu (%d threads) running on %s\n",
            Parallel::get_id(),
            omp_get_num_threads(),
            Parallel::get_host().c_str());
      }
    }
    Parallel::barrier();
  }

  template <class T>
  static void reduce_to_sum(T& t) {
    T t_local = t;
    boost::mpi::all_reduce(Parallel::get_instance().world, t_local, t, std::plus<T>());
  }

  template <class T>
  static void reduce_to_sum_vector(std::vector<T>& t) {
    std::vector<T> t_local = t;
    boost::mpi::all_reduce(Parallel::get_instance().world, t_local, t, vector_plus<T>());
  }
};

#else

// Non-MPI stub for debugging and profiling.
class Parallel {
 public:
  static bool is_master() { return true; }

  static size_t get_id() { return 0; }

  static size_t get_n() { return 1; }

  static std::string get_host() { return "localhost"; }

  static void barrier() {}

  static void print_info() { printf("Running in serial.\n"); }

  template <class T>
  static void reduce_to_sum(T& t) {}
};

#endif  // SERIAL

#endif