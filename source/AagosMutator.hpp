#ifndef AAGOS_MUTATOR_HPP
#define AAGOS_MUTATOR_HPP

#include "AagosOrg.hpp"

#include "emp/math/Distribution.hpp"

namespace aagos {

class AagosMutator {
public:
  using genome_t = AagosOrg::Genome;
  enum class MUTATION_TYPES {
    BIT_FLIPS=0,
    BIT_INSERTIONS,
    BIT_DELETIONS,
    GENE_MOVES,
  };

protected:
  const size_t num_genes;
  const emp::Range<size_t> genome_size_constraints;
  const double prob_gene_moves;
  const double prob_bit_flip;
  const double prob_bit_ins;
  const double prob_bit_del;

  // Used for mutation?
  emp::Binomial gene_moves_binomial;
  emp::vector<emp::Binomial> bit_flips_binomials;
  emp::vector<emp::Binomial> inserts_binomials;
  emp::vector<emp::Binomial> deletes_binomials;

  // Mutation tracking
  std::unordered_map<MUTATION_TYPES, int> last_mutation_tracker;

public:
  AagosMutator(
    size_t n_genes,
    const emp::Range<size_t>& genome_size,
    double p_gene_moves,
    double p_bit_flip,
    double p_bit_ins,
    double p_bit_del
  ) :
    num_genes(n_genes),
    genome_size_constraints(genome_size),
    prob_gene_moves(p_gene_moves),
    prob_bit_flip(p_bit_flip),
    prob_bit_ins(p_bit_ins),
    prob_bit_del(p_bit_del),
    gene_moves_binomial(prob_gene_moves, num_genes)
  {
    for (size_t i = genome_size_constraints.GetLower(); i <= genome_size_constraints.GetUpper(); ++i) {
      bit_flips_binomials.emplace_back(prob_bit_flip, i);
      inserts_binomials.emplace_back(prob_bit_ins, i);
      deletes_binomials.emplace_back(prob_bit_del, i);
    }
  }

  /// Apply gene moves, single-bit substitutions, insertions, and deletions to org's genome.
  size_t ApplyMutations(AagosOrg& org, emp::Random& random) {
    const size_t min_genome_size = genome_size_constraints.GetLower();
    const size_t max_genome_size = genome_size_constraints.GetUpper();
    emp_assert(org.GetNumBits() >= min_genome_size, "Organism's genome exceeds mutator's genome size restrictions.");
    emp_assert(org.GetNumBits() <= max_genome_size, "Organism's genome exceeds mutator's genome size restrictions.");

    auto& genome = org.GetGenome();

    size_t bin_array_offset = org.GetNumBits() - min_genome_size; // offset is num bits - min size of genome
    // Do gene moves
    const size_t num_moves = gene_moves_binomial.PickRandom(random);
    for (size_t m = 0; m < num_moves; ++m) {
      const size_t gene_id = random.GetUInt(0, num_genes);                 // Pick a random gene
      genome.gene_starts[gene_id] = random.GetUInt(genome.bits.GetSize()); // Pick a random new location
    }
    last_mutation_tracker[MUTATION_TYPES::GENE_MOVES] = (int)num_moves;

    // Do bit flips
    const size_t num_flips = bit_flips_binomials[bin_array_offset].PickRandom(random);
    for (size_t m = 0; m < num_flips; ++m) {
      const size_t pos = random.GetUInt(genome.bits.GetSize());
      genome.bits[pos] ^= 1;
    }
    last_mutation_tracker[MUTATION_TYPES::BIT_FLIPS] = (int)num_flips;

    // Do insertions and deletions.
    int num_insert = (int)inserts_binomials[bin_array_offset].PickRandom(random);
    int num_delete = (int)deletes_binomials[bin_array_offset].PickRandom(random);
    const int proj_size = (int)genome.bits.GetSize() + num_insert - num_delete;
    // Check gene size is within range.
    if (proj_size > (int)max_genome_size) {
      num_insert -= (proj_size - (int)max_genome_size);
    } else if (proj_size < (int)min_genome_size) {
      num_delete -= ((int)min_genome_size - proj_size);
    }
    // Assert that we'll be in size limitations.
    emp_assert((int)genome.bits.GetSize() + num_insert - num_delete >= (int)min_genome_size);
    emp_assert((int)genome.bits.GetSize() + num_insert - num_delete <= (int)max_genome_size);
    // Do insertions
    for (int i = 0; i < num_insert; ++i) {
      const size_t pos = random.GetUInt(org.GetNumBits()); // Figure out the position for insertion.
      genome.bits.Resize(genome.bits.GetSize() + 1);       // Increase genome size to make room for insertion.
      emp::BitVector mask(pos, 1);                         // Setup a mask to perserve early bits.
      mask.Resize(genome.bits.GetSize());                     // Align mask size.
      // Now build the new string!
      genome.bits = (mask & genome.bits) | ((genome.bits << 1) & ~mask);
      genome.bits[pos] = random.P(0.5); // Randomize the new bit.
      // Shift any genes that started at pos or later.
      for (auto& x : genome.gene_starts) {
        x += ((size_t)x >= pos);
      }
    }

    // Do deletions
    for (int i = 0; i < num_delete; ++i) {
      size_t pos = random.GetUInt(org.GetNumBits());
      emp::BitVector mask(pos, 1);
      mask.Resize(genome.bits.GetSize());
      genome.bits = (mask & genome.bits) | ((genome.bits >> 1) & ~mask);  // Build the new string!
      genome.bits.Resize(genome.bits.GetSize() - 1);                      // Decrease the size to account for deletion
      // Shift any genes that started at pos or later.
      if (pos == 0) {
        pos = 1; // Adjust position if beginning was deleted (don't want to subtract 1 if gene start = 0)
      }
      for (auto& x : genome.gene_starts) {
        x -= ((size_t)x >= pos);
      }
    }
    last_mutation_tracker[MUTATION_TYPES::BIT_INSERTIONS] = num_insert;
    last_mutation_tracker[MUTATION_TYPES::BIT_DELETIONS] = num_delete;

    // Compute number of mutations, update organism's mutation-related tracking.
    const int num_muts = (int)num_moves + (int)num_flips + num_insert + num_delete;
    emp_assert(num_muts >= 0);
    if (num_muts > 0) {
      org.ResetHistogram();
    }
    return (size_t)num_muts;
  }

  /// Control mutation function. Instead of applying mutations at a per-site rate, apply bit mutations
  /// (i.e., substitutions, insertions, deletions) at a per-gene-per-site rate. This should eliminate
  /// reduced mutational load for compact genetic architectures.
  size_t ApplyMutationsPerGenePerSite(AagosOrg& org, emp::Random& random) {
    const size_t min_genome_size = genome_size_constraints.GetLower();
    const size_t max_genome_size = genome_size_constraints.GetUpper();
    emp_assert(org.GetNumBits() >= min_genome_size, "Organism's genome exceeds mutator's genome size restrictions.");
    emp_assert(org.GetNumBits() <= max_genome_size, "Organism's genome exceeds mutator's genome size restrictions.");

    // This has to be a little more complicated (and less efficient) to take gene occupancy into account when mutating

    genome_t& genome = org.GetGenome();
    // std::cout << "===== Mutation =====" << std::endl;
    // std::cout << "Original (size="<<org.GetGenome().GetNumBits()<<"): ";
    // org.Print(); //
    // std::cout << std::endl;

    // Do gene moves (directly on genome)
    const size_t num_moves = gene_moves_binomial.PickRandom(random);
    for (size_t m = 0; m < num_moves; ++m) {
      const size_t gene_id = random.GetUInt(0, num_genes);                 // Pick a random gene
      // std::cout << "  > moving " << gene_id << std::endl;
      genome.gene_starts[gene_id] = random.GetUInt(genome.bits.GetSize()); // Pick a random new location
    }
    last_mutation_tracker[MUTATION_TYPES::GENE_MOVES] = (int)num_moves;

    // std::cout << "Moves: " << num_moves << std::endl;
    // std::cout << "After moves: "; org.Print();
    // std::cout << std::endl;

    // Build position occupancy map to modify probability of per-site mutations
    const size_t num_genes = genome.num_genes;
    const size_t gene_size = genome.GetGeneSize();
    const size_t genome_size = genome.GetNumBits();

    // Compute gene positions, count occupants per site
    // emp::vector<std::unordered_set> gene_positions(num_genes, {});
    emp::vector<size_t> position_occupants(genome_size, 0);
    for (size_t gene_id = 0; gene_id < num_genes; ++gene_id) {
      const size_t start_pos = genome.gene_starts[gene_id];
      for (size_t pos = start_pos; pos < start_pos + gene_size; ++pos) {
        // gene_positions[gene_id].emplace(pos % genome_size);
        ++position_occupants[pos % genome_size];
      }
    }
    // std::cout << "Position occupant map: " << position_occupants << std::endl;

    // Do bit flips (directly on genome)
    int num_flips = 0;
    for (size_t pos = 0; pos < genome.GetNumBits(); ++pos) {
      const size_t num_occupants = position_occupants[pos];
      bool mutate_bit = random.P(prob_bit_flip); // Apply probability of mutation.
      // for each occupant above the first, apply probability of mutation
      if (num_occupants > 1 && !mutate_bit) {
        for (size_t g = 1; (g < num_occupants) && !mutate_bit; ++g) {
          mutate_bit = random.P(prob_bit_flip);
        }
      }
      if (mutate_bit) {
        genome.bits[pos] ^= 1;
        ++num_flips;
      }
    }
    last_mutation_tracker[MUTATION_TYPES::BIT_FLIPS] = num_flips;
    // std::cout << "Bit flips: " << num_flips << std::endl;
    // std::cout << "After flips: ";
    // org.Print();
    // std::cout << std::endl;

    // Do insertions and deletions (build new genome as we go to maintain accuracy of position_occupants)
    int num_insertions = 0;
    int num_deletions = 0;
    // const size_t start_from = random.GetUInt(genome.GetNumBits()); // Randomize where we start to eliminate
    //                                                                // potential bias in stability of beginning vs end of genome.
    genome_t new_genome(genome);
    new_genome.bits.Resize(genome.GetNumBits() * 2); // Max possible genome growth.
    new_genome.bits.Clear();                    // Reset all bits to 0

    // std::cout << "Making a blank tape (size="<<new_genome.bits.size()<<"): ";
    // new_genome.bits.Print();
    // std::cout << std::endl;
    // std::cout << "blank tape gene starts: " << new_genome.gene_starts << std::endl;

    int write_offset=0; // How much do we offset our writes?
    size_t new_size = genome_size;
    // std::cout << "--insertion/deletions--" << std::endl;
    for (size_t pos = 0; pos < genome_size; ++pos) {
      // std::cout << "pos=" << pos << "; write_offset=" << write_offset << "; new size=" << new_size << std::endl;
      const size_t num_occupants = position_occupants[pos];

      // Do we insert?
      bool do_insertion=random.P(prob_bit_ins); // Apply probability of mutation.
      // for each occupant above the first, apply probability of mutation
      for (size_t g = 1; (g < num_occupants) && !do_insertion; ++g) {
        do_insertion = random.P(prob_bit_ins);
      }
      // Do we delete?
      bool do_deletion=random.P(prob_bit_del);
      // for each occupant above the first, apply probability of mutation
      for (size_t g = 1; (g < num_occupants) && !do_deletion; ++g) {
        do_deletion = random.P(prob_bit_del);
      }

      // std::cout << "  do insertion? " << do_insertion << std::endl;
      // std::cout << "  do deletion? " << do_deletion << std::endl;

      // Do we (1) insert + delete, (2) insert, (3) delete, (4) do nothing
      if (do_insertion && do_deletion) {
        // std::cout << "  > do ins+del" << std::endl;
        // net effect of deletion + insertion is to randomize the bit at this position
        const int write_pos = ((int)pos + write_offset); // should never be negative because pos has to increase 1 for every possible write offset decrement
        emp_assert(write_pos >= 0);
        emp_assert(write_pos < (int)new_genome.bits.size());
        // std::cout << "  > write pos = " << write_pos << std::endl;
        new_genome.bits[write_pos] = random.P(0.5);
        ++num_insertions;
        ++num_deletions;
      } else if (do_insertion && new_size < max_genome_size) {
        // std::cout << "  > do ins" << std::endl;
        const int write_pos = (int)pos + write_offset;
        // std::cout << "  > write pos = " << write_pos << std::endl;
        // insert random bit just before this bit
        new_genome.bits[write_pos] = random.P(0.5);
        // copy original bit
        new_genome.bits[write_pos + 1] = genome.bits.Get(pos);
        // Shift any genes that started at pos or later.
        size_t base_pos=(int)pos + write_offset;
        for (auto & x : new_genome.gene_starts) {
          x += ((size_t)x >= base_pos);
        }
        ++write_offset;
        ++num_insertions;
        ++new_size;
      } else if (do_deletion && new_size > min_genome_size) {
        // std::cout << "  > do del" << std::endl;
        // DON'T copy original bit
        // todo - update gene start positions
        // Shift any genes that started at pos or later.
        size_t base_pos=(int)pos + write_offset;
        if (base_pos==0) base_pos = 1;
        for (auto & x : new_genome.gene_starts) {
          x -= ((size_t)x >= base_pos);
        }
        ++num_deletions;
        --write_offset; // decrement write head offset
        --new_size;
      } else {
        // std::cout << "  > do copy over" << std::endl;
        // copy original bit
        const int write_pos = (int)pos + write_offset;
        // std::cout << "  > write pos = " << write_pos << std::endl;
        new_genome.bits[write_pos] = genome.bits.Get(pos);
      }

      // std::cout << "Updated tape: "; new_genome.bits.Print(); std::cout << std::endl;
      // std::cout << "  updated gene starts: " << new_genome.gene_starts << std::endl;
    }
    emp_assert(new_size >= min_genome_size);
    emp_assert(new_size <= max_genome_size);
    last_mutation_tracker[MUTATION_TYPES::BIT_INSERTIONS] = num_insertions;
    last_mutation_tracker[MUTATION_TYPES::BIT_DELETIONS] = num_deletions;

    // resize new genome
    emp_assert(new_size <= new_genome.bits.size());
    // for (size_t i = 0; i < new_genome.gene_starts.size(); ++i) {
    //   const size_t new_pos = new_genome.gene_starts[i];
    //   const size_t old_pos = genome.gene_starts[i];
    //   emp_assert(new_pos < new_size, genome_size, old_pos, num_insertions, num_deletions, new_pos, new_size);
    // }
    new_genome.bits.Resize(new_size);

    // update organism genome with new genome w/insertions and deletions
    org.GetGenome().bits.Resize(new_genome.bits.size());
    emp_assert(org.GetGenome().gene_starts.size() == new_genome.gene_starts.size());
    for (size_t i = 0; i < org.GetGenome().bits.size(); ++i) {
      org.GetGenome().bits[i] = new_genome.bits.Get(i);
    }
    for (size_t i = 0; i < org.GetGenome().gene_starts.size(); ++i) {
      org.GetGenome().gene_starts[i] = new_genome.gene_starts[i];
    }

    // std::cout << "Final mutated (size="<<org.GetGenome().GetNumBits()<<"): ";
    // org.Print();
    // std::cout << std::endl;
    // std::cout << "  gene starts = " << org.GetGenome().gene_starts << std::endl;

    // Compute number of mutations, update organism's mutation-related tracking.
    const int num_muts = (int)num_moves + (int)num_flips + num_insertions + num_deletions;
    emp_assert(num_muts >= 0);
    if (num_muts > 0) {
      org.ResetHistogram();
    }
    return (size_t)num_muts;
  }

  std::unordered_map<MUTATION_TYPES, int>& GetLastMutations() {
    return last_mutation_tracker;
  }

  void ResetLastMutationTracker() {
    last_mutation_tracker[MUTATION_TYPES::BIT_FLIPS] = 0;
    last_mutation_tracker[MUTATION_TYPES::BIT_INSERTIONS] = 0;
    last_mutation_tracker[MUTATION_TYPES::BIT_DELETIONS] = 0;
    last_mutation_tracker[MUTATION_TYPES::GENE_MOVES] = 0;
  }

};

}

#endif