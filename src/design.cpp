#include "design.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>

#include "utils.h"

namespace SLHD {
Design::Design(const VecInt2D& a, int t) {
  a_ = a;
  n_ = a.size();
  k_ = a[0].size();
  t_ = t;
  m_ = n_ / t;
  assert(m_ * t_ == n_ && n_ > 1);
  if (m_ == 1) {
    m_ = n_;
    t_ = 1;
  }
  InitDis();
  InitVal();
  InitPos();
}

void Design::InitDis() {
  dis_.resize(n_, std::vector<double>(n_));
  for (int i = 0; i < n_; ++i) {
    dis_[i][i] = 0;
    for (int j = 0; j < i; ++j) {
      double d = 0;
      for (int k = 0; k < k_; ++k) {
        d += 1.0 * (a_[i][k] - a_[j][k]) * (a_[i][k] - a_[j][k]);
      }
      dis_[i][j] = dis_[j][i] = std::sqrt(d);
    }
  }
}

void Design::InitVal() {
  sliced_val_.resize(t_, 0);
  total_val_ = 0;
  for (int i = 0; i < n_; ++i) {
    for (int j = 0; j < i; ++j) {
      total_val_ += 1.0 / SLHD::Utils::QuickPow(dis_[i][j], power_);
    }
  }
  total_val_ /= n_ * (n_ - 1) / 2;
  total_val_ = std::pow(total_val_, 1.0 / power_);
  std::vector<double> val_slice(t_, 0);
  for (int s = 0; s < t_; ++s) {
    int l = s * m_;
    int r = (s + 1) * m_;
    for (int i = l; i < r; ++i) {
      for (int j = l; j < i; ++j) {
        sliced_val_[s] += 1.0 / SLHD::Utils::QuickPow(dis_[i][j], power_);
      }
    }
    sliced_val_[s] /= m_ * (m_ - 1) / 2;
    sliced_val_[s] = std::pow(sliced_val_[s], 1.0 / power_);
  }
  double sub_val = std::accumulate(sliced_val_.begin(), sliced_val_.end(), 0.0);
  val_ = 0.5 * (total_val_ + sub_val / t_);
}

void Design::InitPos() {
  pos_.resize(k_, VecInt2D(m_, std::vector<int>(t_)));
  for (int i = 0; i < n_; ++i) {
    for (int j = 0; j < k_; ++j) {
      int s = i / m_;
      int l = (a_[i][j] - 1) / t_;
      pos_[j][l][s] = i;
    }
  }
}

void Design::SwapInCol(int col, int r1, int r2) {
  auto Update = [this, col](int r1, int r2, int delta, double& val, bool update_dis = false) {
    auto GetNewDis = [this, col](int r1, int r2, int delta) -> double {
      double tmp_dis = dis_[r1][r2] * dis_[r1][r2];
      tmp_dis -= 1.0 * (a_[r2][col] - a_[r1][col]) * (a_[r2][col] - a_[r1][col]);
      tmp_dis += 1.0 * (a_[r2][col] - a_[r1][col] + delta) * (a_[r2][col] - a_[r1][col] + delta);
      return std::sqrt(tmp_dis);
    };
    val -= 1.0 / Utils::QuickPow(dis_[r1][r2], power_);
    double new_dis = GetNewDis(r1, r2, delta);
    val += 1.0 / Utils::QuickPow(new_dis, power_);
    if (update_dis) dis_[r1][r2] = new_dis;
  };
  // Update sliced value
  if (r1 / m_ == r2 / m_) {
    int slice_id = r1 / m_;
    double sliced_tmp_val = Utils::QuickPow(sliced_val_[slice_id], power_);
    sliced_tmp_val *= m_ * (m_ - 1) / 2;
    int delta = a_[r2][col] - a_[r1][col];
    for (int i = slice_id * m_; i < (slice_id + 1) * m_; ++i) {
      if (i == r1 || i == r2) continue;
      i > r1 ? Update(i, r1, delta, sliced_tmp_val) : Update(r1, i, -delta, sliced_tmp_val);
      i > r2 ? Update(i, r2, -delta, sliced_tmp_val) : Update(r2, i, delta, sliced_tmp_val);
    }
    sliced_tmp_val /= m_ * (m_ - 1) / 2;
    sliced_tmp_val = std::pow(sliced_tmp_val, 1.0 / power_);
    val_ += (sliced_tmp_val - sliced_val_[slice_id]) / t_ / 2;
    sliced_val_[slice_id] = sliced_tmp_val;
  } else {
    int delta = a_[r2][col] - a_[r1][col];
    int s1 = r1 / m_;
    double s1_tmp_val = Utils::QuickPow(sliced_val_[s1], power_);
    s1_tmp_val *= m_ * (m_ - 1) / 2;
    for (int i = s1 * m_; i < (s1 + 1) * m_; ++i) {
      if (i == r1) continue;
      i > r1 ? Update(i, r1, delta, s1_tmp_val) : Update(r1, i, -delta, s1_tmp_val);
    }
    s1_tmp_val /= m_ * (m_ - 1) / 2;
    s1_tmp_val = std::pow(s1_tmp_val, 1.0 / power_);
    val_ += (s1_tmp_val - sliced_val_[s1]) / t_ / 2;
    sliced_val_[s1] = s1_tmp_val;

    int s2 = r2 / m_;
    double s2_tmp_val = Utils::QuickPow(sliced_val_[s2], power_);
    s2_tmp_val *= m_ * (m_ - 1) / 2;
    for (int i = s2 * m_; i < (s2 + 1) * m_; ++i) {
      if (i == r2) continue;
      i > r2 ? Update(i, r2, -delta, s2_tmp_val) : Update(r2, i, delta, s2_tmp_val);
    }
    s2_tmp_val /= m_ * (m_ - 1) / 2;
    s2_tmp_val = std::pow(s2_tmp_val, 1.0 / power_);
    val_ += (s2_tmp_val - sliced_val_[s2]) / t_ / 2;
    sliced_val_[s2] = s2_tmp_val;
  }
  // Update total value
  double total_tmp_val = Utils::QuickPow(total_val_, power_);
  total_tmp_val *= n_ * (n_ - 1) / 2;
  int delta = a_[r2][col] - a_[r1][col];
  for (int i = 0; i < n_; ++i) {
    if (i == r1 || i == r2) continue;
    i > r1 ? Update(i, r1, delta, total_tmp_val, true) : Update(r1, i, -delta, total_tmp_val, true);
    i > r2 ? Update(i, r2, -delta, total_tmp_val, true) : Update(r2, i, delta, total_tmp_val, true);
  }
  total_tmp_val /= n_ * (n_ - 1) / 2;
  total_tmp_val = std::pow(total_tmp_val, 1.0 / power_);
  val_ += (total_tmp_val - total_val_) / 2;
  total_val_ = total_tmp_val;

  // Update position
  if (r1 / m_ != r2 / m_) {
    assert((a_[r1][col] - 1) / t_ == (a_[r2][col] - 1) / t_); // same level swap between two slices
    pos_[col][(a_[r1][col] - 1) / t_][r2 / m_] = r2;
    pos_[col][(a_[r2][col] - 1) / t_][r1 / m_] = r1;
  }

  // Update design matrix
  std::swap(a_[r1][col], a_[r2][col]);
}

void Design::Display() const {
  std::cout << m_ << " " << k_ << " " << t_ << "\n";
  int w = std::floor(std::log10(n_)) + 1;
  for (int i = 0; i < n_; ++i) {
    for (int j = 0; j < k_; ++j) {
      std::cout << std::setw(w) << a_[i][j] << " \n"[j == k_ - 1];
    }
  }
  std::cout << "Criterion value: " << val_ << std::endl;
}
} // namespace SLHD