#pragma once

#include <string>
#include <vector>

namespace SLHD {
class Design {
 public:
  using VecInt2D = std::vector<std::vector<int>>;
  using VecInt3D = std::vector<VecInt2D>;
  using VecDouble2D = std::vector<std::vector<double>>;
  Design(const VecInt2D& a, int t = 1);
  inline int GetN() const { return n_; }
  inline int GetM() const { return m_; }
  inline int GetK() const { return k_; }
  inline int GetT() const { return t_; }
  inline const VecInt2D& GetA() const { return a_; }
  inline const VecDouble2D& GetDis() const { return dis_; }
  inline double GetVal() const { return val_; }
  inline const VecInt3D& GetPos() const { return pos_; }
  inline void SetPower(int power) { power_ = power; }
  void SwapInCol(int col, int r1, int r2);
  void Display() const;
 private:
  void InitDis();
  void InitVal();
  void InitPos();
 private:
  int power_ = 15;   // Parameter p in phip criterion, 15 by default
  int m_;            // Number of runs of slices
  int k_;            // Number of (quantitative) factors
  int t_;            // Number of slices
  int n_;            // Number of total runs
  VecInt2D a_;       // dim (n, k)
  VecDouble2D dis_;  // Euclidean distance matrix, dim (n, n)
  double total_val_; // Phip value of LHD
  std::vector<double> sliced_val_; // Phip value of each slices, dim (t,)
  double val_;       // Phip value of SLHD
  VecInt3D pos_;     // Row position of each level, dim (k, m, t)
};
} // namespace SLHD