#include "manager.h"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <random>

namespace SLHD {
Manager::Manager(int m, int k, int t)
    : m_(m),
      k_(k),
      t_(t),
      n_(m * t),
      solver_(m, k, t) {
  ranges_.resize(k_, {0, 1});
}

void Manager::Display(const std::vector<std::vector<double>>& design) {
  std::cout << std::fixed << std::setprecision(2);
  double max_num = 0;
  for (int i = 0; i < design.size(); ++i) {
    max_num = std::max(max_num, *std::max_element(design[i].begin(), design[i].end()));
  }
  int w = std::floor(std::log10(max_num)) + 1 + 4;
  for (int i = 0; i < design.size(); ++i) {
    for (int j = 0; j < design[i].size(); ++j) {
      std::cout << std::setw(w) << design[i][j] << " \n"[j + 1 == design[i].size()];
    }
  }
  std::cout.flush();
}

void Manager::SetQuantitativeFactorRange(const std::vector<std::pair<double, double>>& ranges) {
  assert(ranges.size() == k_);
  ranges_ = ranges;
}

void Manager::SetQualitativeFactorLevel(const std::vector<std::vector<double>>& levels) {
  assert(levels.size() == t_);
  for (int i = 1; i < t_; ++i) {
    assert(levels[i].size() == levels[i - 1].size());
  }
  levels_ = levels;
}

std::vector<std::vector<double>> Manager::GetDesign() {
  std::mt19937 rng(seed_);
  solver_.SetSeed(seed_);
  std::vector<std::vector<int>> ori_design = solver_.GetOptimalSLHD();
  int factor_num = k_ + (levels_.empty() ? 0 : levels_[0].size());
  std::vector<std::vector<double>> design(n_, std::vector<double>(factor_num));
  std::uniform_real_distribution<double> uniform_dis(std::nextafter(0.0, 1.0), 1.0);
  for (int i = 0; i < n_; ++i) {
    // fill quantitative factors
    for (int j = 0; j < k_; ++j) {
      double e = uniform_dis(rng);
      design[i][j] = (ori_design[i][j] - e) / n_;
      design[i][j] = ranges_[j].first +
        (ranges_[j].second - ranges_[j].first) * design[i][j];
    }
    // fill qualitative factors
    for (int j = k_; j < factor_num; ++j) {
      int slice_id = i / m_;
      design[i][j] = levels_[slice_id][j - k_];
    }
  }
  return design;
}
} // namespace SLHD