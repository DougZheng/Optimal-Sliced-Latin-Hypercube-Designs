#pragma once

#include <iostream>
#include <random>

#include "design.h"

// Optimal Sliced Latin Hypercube Designs
// Ref: https://www.asc.ohio-state.edu/statistics/comp_exp/jour.club/optimal_sliced_lhd_ba2015.pdf

namespace SLHD {
class SLHD {
 public:
  SLHD(int m, int k, int t = 1);
  static void Display(const std::vector<std::vector<int>>& design, int t = 1);
  inline void SetIMax(int imax) { imax_ = imax; }
  inline void SetIterCnt(int iter_cnt) { iter_cnt_ = iter_cnt; }
  inline void SetInitTemp(int temp_init) { temp_init_ = temp_init; }
  inline void SetTempRate(int temp_rate) { temp_rate_ = temp_rate; }
  inline void SetSeed(int seed) { rng_.seed(seed); }
  inline void SetLogFreq(int log_freq) { log_freq_ = log_freq; }
  inline void SetLogOp(bool log_op) { log_op_ = log_op; }
  std::vector<std::vector<int>> GetRandomSLHD();
  std::vector<std::vector<int>> GetOptimalSLHD();
 private:
  inline int Rand2(int n, int x) {
    int y = rng_() % n;
    while (x == y) y = rng_() % n;
    return y;
  }
  inline void Log(int cnt, const Design& design) const {
    if (!log_op_ || cnt % log_freq_ != 0) return;
    std::cerr << "Iteration: " << cnt << ", " <<
      "Val: " << design.GetVal() << std::endl;
  }
  void InitDefaultParam();
  Design InnerGetOptimalSLHD(Design design, double p0, int max_iter_cnt);
 private: 
  int m_;             // Number of runs of slices
  int k_;             // Number of (quantitative) factors
  int t_;             // Number of slices
  int n_;             // Number of total runs
  int imax_;          // Maximum number of tries without acceptance
  int iter_cnt_;      // Maximum number of iterations
  double temp_init_;  // Initial temperature
  double temp_rate_;  // Annealing rate
  std::mt19937 rng_;  // Random number engine
  int log_freq_;      // Log frequency
  bool log_op_;       // Log option
};
} // namespace SLHD