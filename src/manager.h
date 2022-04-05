#pragma once

#include <vector>
#include <utility>

#include "SLHD.h"

namespace SLHD {
class Manager {
 public:
  Manager(int m, int k, int t = 1);
  // Fix seed to fix the resulting design
  inline void SetSeed(int seed) { seed_ = seed; }
  // Get the reference of SLHD solver to change parameters of the algorithm
  // See the definition of class SLHD
  inline SLHD& GetSolverRef() { return solver_; }
  static void Display(const std::vector<std::vector<double>>& design);
  // Default range of each quantitative factors is [0, 1]
  void SetQuantitativeFactorRange(const std::vector<std::pair<double, double>>& ranges);
  // No qualitative factors by default
  void SetQualitativeFactorLevel(const std::vector<std::vector<double>>& levels);
  std::vector<std::vector<double>> GetDesign();
 private:
  int m_;       // Number of runs of slices
  int k_;       // Number of quantitative factors
  int t_;       // Number of slices (number of level combinations of qualitative factors)
  int n_;       // Number of total runs
  int seed_;    // Seed of random number engine
  SLHD solver_; // SLHD solver

  // Ranges of quantitative factors, dim (k,)
  std::vector<std::pair<double, double>> ranges_;  
  // Level combinations of qualitative factors, 
  // dim (t, c) where c is the number of qualitative factors
  std::vector<std::vector<double>> levels_;
};
};