#ifndef TIME_H_
#define TIME_H_

#include "parallel.h"
#include "std.h"

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_RESET "\x1b[0m"

class Time {
 public:
  // To be called at the beginning of the program.
  static void init() {
    Parallel::barrier();
    if (Parallel::is_master()) {
      const auto& now = std::chrono::high_resolution_clock::now();
      Time::get_instance().init_time = now;
      Time::get_instance().prev_time = now;
    }
  }

  static void start(const std::string& event) {
    Parallel::barrier();
    if (Parallel::is_master()) {
      const auto& now = std::chrono::high_resolution_clock::now();
      const auto& init_time = Time::get_instance().init_time;
      auto& start_times = Time::get_instance().start_times;
      start_times.push_back(std::make_pair(event, now));
      printf("\n");
      printf("START ");
      for (size_t i = 0; i < start_times.size() - 1; i++) {
        printf("%s >> ", start_times[i].first.c_str());
      }
      printf(
          "%s " ANSI_COLOR_YELLOW "[%.3f]" ANSI_COLOR_RESET "\n",
          event.c_str(),
          get_duration(init_time, now));
      Time::get_instance().prev_time = now;
    }
    Parallel::barrier();
  }

  static void end() {
    Parallel::barrier();
    if (Parallel::is_master()) {
      const auto& now = std::chrono::high_resolution_clock::now();
      const auto& init_time = Time::get_instance().init_time;
      auto& start_times = Time::get_instance().start_times;
      printf("==END ");
      for (size_t i = 0; i < start_times.size() - 1; i++) {
        printf("%s >> ", start_times[i].first.c_str());
      }
      const auto& event_name = start_times.back().first;
      const auto& event_start_time = start_times.back().second;
      printf(
          "%s " ANSI_COLOR_YELLOW "[%.3f/%.3f]" ANSI_COLOR_RESET "\n",
          event_name.c_str(),
          get_duration(event_start_time, now),
          get_duration(init_time, now));
      start_times.pop_back();
      Time::get_instance().prev_time = now;
    }
    Parallel::barrier();
  }

  static void checkpoint(const std::string& msg) {
    Parallel::barrier();
    if (Parallel::is_master()) {
      const auto& now = std::chrono::high_resolution_clock::now();
      const auto& init_time = Time::get_instance().init_time;
      auto& prev_time = Time::get_instance().prev_time;
      auto& start_times = Time::get_instance().start_times;
      const auto& event_start_time = start_times.back().second;
      printf(
          "--CHK %s " ANSI_COLOR_YELLOW "[%.3f/%.3f/%.3f]" ANSI_COLOR_RESET "\n",
          msg.c_str(),
          get_duration(prev_time, now),
          get_duration(event_start_time, now),
          get_duration(init_time, now));
      prev_time = now;
    }
    Parallel::barrier();
  }

 private:
  std::chrono::high_resolution_clock::time_point init_time;
  std::chrono::high_resolution_clock::time_point prev_time;
  std::vector<std::pair<std::string, std::chrono::high_resolution_clock::time_point>> start_times;

  // Singleton pattern.
  static Time& get_instance() {
    static Time instance;
    return instance;
  }

  static double get_duration(
      const std::chrono::high_resolution_clock::time_point start,
      const std::chrono::high_resolution_clock::time_point end) {
    return (std::chrono::duration_cast<std::chrono::duration<double>>(end - start)).count();
  }
};

#endif