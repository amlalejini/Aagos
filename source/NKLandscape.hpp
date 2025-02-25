#ifndef AAGOS_NKLANDSCAPE_H
#define AAGOS_NKLANDSCAPE_H

#include "emp/base/vector.hpp"
#include "emp/base/assert.hpp"
#include "emp/math/Random.hpp"
#include "emp/math/math.hpp"
#include "emp/bits/BitVector.hpp"

/*
  NOTE: This class is adapted from the NKLandscape class in the Empirical library
    - https://github.com/devosoft/Empirical/blob/master/include/emp/Evolve/NK.hpp
*/

namespace aagos {

/// An NK Landscape is a popular tool for studying theoretical questions about evolutionary
/// dynamics. It is a randomly generated fitness landscape on which bitstrings can evolve.
/// NK Landscapes have two parameters: N (the length of the bitstrings) and K (epistasis).
/// Since you have control over the amount of epistasis, NK Landscapes are often called
/// "tunably rugged" - a useful feature, since the ruggedness of the fitness landscape is thought
/// to be important to many evolutionary dynamics. For each possible value that a site and its
/// K neighbors to the right can have, a random fitness contribution is chosen.
/// These contributions are summed across the bitstring. So when K = 0, each site has a single
/// optimal value, resulting in a single smooth fitness peak.
///
/// For more information, see Kauffman and Levin,
/// 1987 (Towards a general theory of adaptive walks on rugged landscapes).
///
/// This object handles generating and maintaining an NK fitness landscape.
/// Note: Overly large Ns and Ks currently trigger a seg-fault, caused by trying to build a table
/// that is larger than will fit in memory. If you are using small values for N and K,
/// you can get better performance by using an NKLandscapeConst instead.

class NKLandscape {
  private:
    size_t N;             ///< The number of bits in each genome.
    size_t K;             ///< The number of OTHER bits with which each bit is epistatic.
    size_t state_count;   ///< The total number of states associated with each bit table.
    size_t total_count;   ///< The total number of states in the entire landscape space.
    emp::vector< emp::vector<double> > landscape;  ///< The actual values in the landscape.

  public:
    NKLandscape() : N(0), K(0), state_count(0), total_count(0), landscape() { ; }
    NKLandscape(const NKLandscape &) = default;
    NKLandscape(NKLandscape &&) = default;

    /// N is the length of bitstrings in your population, K is the number of neighboring sites
    /// the affect the fitness contribution of each site (i.e. epistasis or ruggedness), random
    /// is the random number generator to use to generate this landscape.
    NKLandscape(size_t _N, size_t _K, emp::Random & random)
      : N(_N), K(_K)
      , state_count(emp::IntPow<size_t>(2,K+1))
      , total_count(N * state_count)
      , landscape(N)
    {
      Reset(random);
    }
    ~NKLandscape() { ; }
    NKLandscape & operator=(const NKLandscape &) = delete;
    NKLandscape & operator=(NKLandscape &&) = default;

    /// Randomize the landscape without changing the landscape size.
    void Reset(emp::Random & random) {
      emp_assert(K < 32, K);
      emp_assert(K < N, K, N);

      // Build new landscape.
      for ( auto & ltable : landscape) {
        ltable.resize(state_count);
        for (double & pos : ltable) {
          pos = random.GetDouble();
        }
      }
    }

    /// Configure for new values of N and K.
    void Config(size_t _N, size_t _K, emp::Random & random) {
      // Save new values.
      N = _N;  K = _K;
      state_count = emp::IntPow<size_t>(2,K+1);
      total_count = N * state_count;
      landscape.resize(N);
      Reset(random);
    }

    /// Returns N
    size_t GetN() const { return N; }
    /// Returns K
    size_t GetK() const { return K; }
    /// Get the number of posssible states for a given site
    size_t GetStateCount() const { return state_count; }
    /// Get the total number of states possible in the landscape
    /// (i.e. the number of different fitness contributions in the table)
    size_t GetTotalCount() const { return total_count; }

    const emp::vector<emp::vector<double>>& GetLandscape() const { return landscape; }

    /// Get the fitness contribution of position [n] when it (and its K neighbors) have the value
    /// [state]
    double GetFitness(size_t n, size_t state) const {
      emp_assert(state < state_count, state, state_count);
      return landscape[n][state];
    }

    /// Get the fitness of a whole  bitstring
    double GetFitness( std::vector<size_t> states ) const {
      emp_assert(states.size() == N);
      double total = landscape[0][states[0]];
      for (size_t i = 1; i < N; i++) total += GetFitness(i,states[i]);
      return total;
    }

    /// Get the fitness of a whole bitstring (pass by value so can be modified.)
    double GetFitness(emp::BitVector genome) const {
      emp_assert(genome.GetSize() == N, genome.GetSize(), N);

      // Use a double-length genome to easily handle wrap-around.
      genome.Resize(N*2);
      genome |= (genome << N);

      double total = 0.0;
      size_t mask = emp::MaskLow<size_t>(K+1);
      for (size_t i = 0; i < N; i++) {
        const size_t cur_val = (genome >> i).GetUInt(0) & mask;
        const double cur_fit = GetFitness(i, cur_val);
        total += cur_fit;
      }
      return total;
    }

    /// Get the fitness of a site in a bitstring
    // (pass by value so can be modified.)
    double GetSiteFitness(size_t n, emp::BitVector genome) const {
      emp_assert(genome.GetSize() == N, genome.GetSize(), N);

      // Use a double-length genome to easily handle wrap-around.
      genome.Resize(N*2);
      genome |= (genome << N);

      size_t mask = emp::MaskLow<size_t>(K+1);
      const size_t cur_val = (genome >> n).GetUInt(0) & mask;
      return GetFitness(n, cur_val);
    }

    /// Get the fitness of a whole bitstring (pass by value so can be modified.)
    emp::vector<double> GetFitnesses(emp::BitVector genome) const {
      // Use a double-length genome to easily handle wrap-around.
      genome.Resize(N*2);
      genome |= (genome << N);
      emp::vector<double> fits;

      size_t mask = emp::MaskLow<size_t>(K+1);
      for (size_t i = 0; i < N; i++) {
        const size_t cur_val = (genome >> i).GetUInt(0) & mask;
        const double cur_fit = GetFitness(i, cur_val);
        fits.push_back(cur_fit);
      }
      return fits;
    }

    void SetState(size_t n, size_t state, double in_fit) { landscape[n][state] = in_fit; }

    void RandomizeStates(emp::Random & random, size_t num_states=1) {
      for (size_t i = 0; i < num_states; i++) {
        SetState(random.GetUInt(N), random.GetUInt(state_count), random.GetDouble());
      }
    }

};

} // End aagos namespace

#endif