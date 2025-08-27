#ifndef GRADIENT_FITNESS_MODEL_HPP
#define GRADIENT_FITNESS_MODEL_HPP

#include "emp/base/vector.hpp"
#include "emp/bits/Bits.hpp"
#include "emp/math/random_utils.hpp"
#include "emp/tools/string_utils.hpp"

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

namespace aagos {

/// Fitness model for gradient fitness evaluation
struct GradientFitnessModel {
  size_t num_genes;
  size_t gene_size;
  emp::vector<emp::BitVector> targets;

  GradientFitnessModel(emp::Random & rand, size_t n_genes, size_t g_size)
    : num_genes(n_genes), gene_size(g_size)
  {
    for (size_t i = 0; i < num_genes; ++i) {
      targets.emplace_back(emp::RandomBitVector(rand, gene_size));
      emp_assert(targets.back().GetSize() == gene_size);
    }
    emp_assert(targets.size() == num_genes);
  }

  const emp::BitVector & GetTarget(size_t id) const { return targets[id]; }

  /// Mutate a number of target bits equal to bit cnt.
  void RandomizeTargetBits(emp::Random & rand, size_t bit_cnt) {
    for (size_t i = 0; i < bit_cnt; ++i) {
      // Select a random target sequence.
      const size_t target_id = rand.GetUInt(targets.size());
      emp::BitVector & target = targets[target_id];
      const size_t target_pos = rand.GetUInt(target.GetSize());
      target.Set(target_pos, !target.Get(target_pos));
    }
  }

  /// Randomize a number of targets equal to cnt.
  void RandomizeTargets(emp::Random & rand, size_t cnt) {
    // Change a number of targets (= cnt).
    emp::vector<size_t> target_ids(targets.size());
    std::iota(target_ids.begin(), target_ids.end(), 0);
    emp::Shuffle(rand, target_ids);
    cnt = emp::Min(cnt, target_ids.size()); // No need to randomize more targets than exist.
    for (size_t i = 0; i < cnt; ++i) {
      // Select a random target sequence.
      emp_assert(i < target_ids.size());
      const size_t target_id = target_ids[i];
      emp::BitVector & target = targets[target_id];
      emp::RandomizeBitVector(target, rand);
    }
  }

  void PrintTargets(std::ostream & out=std::cout) {
    // lazily outsource to emp::to_string
    out << emp::to_string(targets);
  }

  /// Load targets from file (specified by given path)
  /// First non-commented line of file should give what emp::to_string would output.
  ///  - [ bits bits bits bits ]
  ///  - Loaded environment must be consistent with num_genes and gene_size
  bool LoadTargets(const std::string & path) {
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
      if (cur_line == emp::empty_string()) continue; // Skip empty line.
      else if (cur_line[0] == '#') continue; // Treat '#' as a commented line
      else if (cur_line[0] == '[') {
        // Attempt to read environments state.
        emp::remove_chars(cur_line, "[]"); // Remove brackets
        emp::left_justify(cur_line);       // Remove leading whitespace
        emp::right_justify(cur_line);      // Remove trailing whitespace
        line_components.clear();
        emp::slice(cur_line, line_components, ' ');
        // Each slice should be a gene target, check to make sure number of targets is correct.
        if (line_components.size() != num_genes) return false;

        for (size_t i = 0; i < line_components.size(); ++i) {
          emp_assert(i < targets.size());
          auto & target_str = line_components[i];
          // Target string should be correct size.
          if (target_str.size() != gene_size) return false;
          for (size_t bit = 0; bit < target_str.size(); ++bit) {
            emp_assert(bit < targets[i].GetSize());
            // if (component == '1')
            targets[i].Set(targets[i].GetSize() - bit - 1, target_str[bit] == '1');
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