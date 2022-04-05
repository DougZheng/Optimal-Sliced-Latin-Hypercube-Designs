#include "SLHD.h"

#include <iostream>
#include <string>

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage: ${exe} ${m} ${k} [${t} ${iter_cnt}]" << std::endl;
    return 1;
  }
  int m = std::stoi(argv[1]);
  int k = std::stoi(argv[2]);
  int t = argc > 3 ? std::stoi(argv[3]) : 1;
  SLHD::SLHD solver(m, k, t);
  if (argc > 4) {
    int iter_cnt = std::stoi(argv[4]);
    solver.SetIterCnt(iter_cnt);
  }
  solver.SetLogOp(true);
  auto design = solver.GetOptimalSLHD();
  solver.Display(design, t);
  return 0;
}