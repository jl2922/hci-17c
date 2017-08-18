#ifndef CONFIG_H_
#define CONFIG_H_

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include "std.h"

class Config {
 private:
  boost::property_tree::ptree config_tree;

  Config(){};  // Prevent instantiation.

  // Singleton pattern.
  static Config& get_instance() {
    static Config config;
    return config;
  }

  static boost::property_tree::ptree& get_config_tree() {
    return Config::get_instance().config_tree;
  }

 public:
  static void load(const std::string& filename) {
    boost::property_tree::read_json(filename, Config::get_config_tree());
  }

  template <class T>
  static T get(const std::string& property) {
    return Config::get_config_tree().get<T>(property);
  }

  template <class T>
  static T get(const std::string& property, const T& default_value) {
    return Config::get_config_tree().get<T>(property, default_value);
  }

  template <class T>
  static std::vector<T> get_array(const std::string& property) {
    std::vector<T> res;
    for (auto& item : Config::get_config_tree().get_child(property)) {
      res.push_back(item.second.get_value<T>());
    }
    return res;
  }

  static void print() {
    printf("Configuration:\n");
    boost::property_tree::json_parser::write_json(std::cout, Config::get_config_tree());
    printf("\n");
  }
};

#endif