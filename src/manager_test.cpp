#include "manager.h"

#include <ctime>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <utility>
#include <random>

std::vector<std::pair<double, double>> GenRanges(int k) {
  std::mt19937 rng(std::time(nullptr));
  std::vector<std::pair<double, double>> ranges(k);
  std::cout << "Quantitative factor ranges:\n";
  for (int i = 0; i < k; ++i) {
    double l = rng() % 100;
    double r = rng() % 100;
    if (l > r) std::swap(l, r);
    ranges[i] = {l, r};
    std::cout << "Factor " << i << ", range [" <<
      std::setw(2) << l << ", " << std::setw(2) << r << "]\n";
  }
  std::cout.flush();
  return ranges;
};

std::vector<std::vector<double>> GenLevels(int t) {
  std::mt19937 rng(std::time(nullptr));
  int factor_num = rng() % 3 + 1;
  std::vector<std::vector<double>> levels(t, std::vector<double>(factor_num));
  std::cout << "Qualitative factor combination levels:\n";
  for (int i = 0; i < t; ++i) {
    std::cout << "Level " << i << ":";
    for (int j = 0; j < factor_num; ++j) {
      levels[i][j] = rng() % 100;
      std::cout << std::setw(3) << levels[i][j] << ",\n"[j == factor_num - 1];
    }
  }
  std::cout.flush();
  return levels;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage: ${exe} ${m} ${k} [${t} ${iter_cnt}]" << std::endl;
    return 1;
  }
  int m = std::stoi(argv[1]);
  int k = std::stoi(argv[2]);
  int t = argc > 3 ? std::stoi(argv[3]) : 1;
  SLHD::Manager slhd_manager(m, k, t);
  slhd_manager.SetSeed(20220405); // Optional
  slhd_manager.SetQuantitativeFactorRange(GenRanges(k));
  slhd_manager.SetQualitativeFactorLevel(GenLevels(t));
  {
    // Optional, setting parameters of the algorithm
    auto& slhd_solver = slhd_manager.GetSolverRef();
    if (argc > 4) {
      int iter_cnt = std::stoi(argv[4]);
      slhd_solver.SetIterCnt(iter_cnt);
    }
    slhd_solver.SetLogOp(true);
  }
  auto design = slhd_manager.GetDesign();
  slhd_manager.Display(design);
  return 0;
}