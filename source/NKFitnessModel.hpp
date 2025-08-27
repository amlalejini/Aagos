#ifndef NK_FITNESS_MODEL_HPP
#define NK_FITNESS_MODEL_HPP

#include "NKLandscape.hpp"

#include "emp/base/vector.hpp"
#include "emp/bits/Bits.hpp"
#include "emp/math/random_utils.hpp"
#include "emp/tools/string_utils.hpp"

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

namespace aagos {

/// Fitness model for NK fitness evaluation
struct NKFitnessModel {
  size_t num_genes;
  size_t gene_size;
  aagos::NKLandscape landscape;

  NKFitnessModel(emp::Random& rand, size_t n_genes, size_t g_size)
    : num_genes(n_genes), gene_size(g_size)
  {
    landscape.Config(num_genes, gene_size - 1, rand);
  }

  aagos::NKLandscape& GetLandscape() { return landscape; }

  void RandomizeLandscapeBits(emp::Random& rand, size_t cnt) {
    landscape.RandomizeStates(rand, cnt);
  }

  void PrintLandscape(std::ostream& out=std::cout) {
    out << emp::to_string(landscape.GetLandscape());
  }

  /// Load NK landscape from file (specified by given path)
  /// First non-commented line of file should give what emp::to_string would output.
  ///  - [ [ fitness fitness fitness ... ] [ fitness fitness fitness ... ] ... ]
  ///  - Loaded landscape must be consistent with num_genes and gene_size
  bool LoadLandscape(const std::string & path) {
    std::ifstream environment_fstream(path);
    if (!environment_fstream.is_open()) {
      std::cout << "Failed to open environment file (" << path << "). Exiting..." << std::endl;
      exit(-1);
    }
    std::string cur_line;
    emp::vector<std::string> line_components;
    bool success = false;
    while (!environment_fstream.eof()) {
      std::getline(environment_fstream, cur_line);
      emp::left_justify(cur_line); // Remove leading whitespace
      if (cur_line == emp::empty_string()) continue; // Skip empty lines
      else if (cur_line[0] == '#') continue; // Treat '#' as a commented line
      else if (cur_line[0] == '[') {
        // Attempt to read landscape state.
        emp::remove_chars(cur_line, "["); // Remove opening brackets. Using closing brackets to delimit gene-associated states.
        emp::left_justify(cur_line);
        emp::right_justify(cur_line);
        line_components.clear();
        emp::slice(cur_line, line_components, ']');
        // Last line component is blank space.
        line_components.pop_back();
        if (landscape.GetN() != line_components.size()) return false;
        for (size_t n = 0; n < line_components.size(); ++n) {
          std::string & values_str = line_components[n];
          emp::vector<std::string> values;
          emp::left_justify(values_str);
          emp::right_justify(values_str);
          emp::slice(values_str, values, ' ');
          if (landscape.GetStateCount() != values.size()) return false;
          for (size_t state = 0; state < values.size(); ++state) {
            double value = emp::from_string<double>(values[state]);
            landscape.SetState(n, state, value);
          }
        }
        success = true;
        break;
      }
    }
    return success;
  }
};

}

#endif