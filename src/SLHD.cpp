#include "SLHD.h"

#include <cmath>
#include <algorithm>

namespace SLHD {
SLHD::SLHD(int m, int k, int t) : m_(m), k_(k), t_(t), n_(m * t) {
  // Equivalent conversion to get maximinLHD
  if (m_ == 1) {
    m_ = n_;
    t_ = 1;
  }
  InitDefaultParam();
}

void SLHD::InitDefaultParam() {
  imax_ = 100;
  iter_cnt_ = 100000;
  temp_rate_ = 0.05;
  rng_.seed(0);
  double avg_d2 = 1.0 * k_ * n_ * (n_ + 1) / 6;
  double delta = (1.0 / std::sqrt(avg_d2 - k_) - 1.0 / std::sqrt(avg_d2));
  temp_init_ = -delta / std::log(0.99);
  log_freq_ = 10000;
  log_op_ = false;
}

void SLHD::Display(const std::vector<std::vector<int>>& design, int t) {
  Design(design, t).Display();
}

std::vector<std::vector<int>> SLHD::GetRandomSLHD() {
  std::vector<std::vector<int>> a(n_, std::vector<int>(k_));
  std::vector<int> perm(t_);
  std::iota(perm.begin(), perm.end(), 1);
  for (int i = 0; i < k_; ++i) {
    for (int j = 0; j < m_; ++j) {
      std::shuffle(perm.begin(), perm.end(), rng_);
      for (int k = 0; k < t_; ++k) {
        a[k * m_ + j][i] = j * t_ + perm[k];
      }
    }
  }
  std::vector<int> nums(m_);
  for (int i = 0; i < t_; ++i) {
    for (int j = 0; j < k_; ++j) {
      for (int k = 0; k < m_; ++k) {
        nums[k] = a[i * m_ + k][j];
      }
      std::shuffle(nums.begin(), nums.end(), rng_);
      for (int k = 0; k < m_; ++k) {
        a[i * m_ + k][j] = nums[k];
      }
    }
  }
  return a;
}

std::vector<std::vector<int>> SLHD::GetOptimalSLHD() {
  if (n_ == 1) {
    return std::vector<std::vector<int>>(1, std::vector<int>(k_, 1));
  }
  auto design = Design(GetRandomSLHD(), t_);
  int state1_iter_cnt = t_ == 1 ? iter_cnt_ : static_cast<int>(std::ceil(iter_cnt_ * 0.75));
  design = InnerGetOptimalSLHD(design, 1.0, state1_iter_cnt);
  if (t_ > 1) {
    int state2_iter_cnt = iter_cnt_ - state1_iter_cnt;
    design = InnerGetOptimalSLHD(design, 0.0, state2_iter_cnt);
  }
  return design.GetA();
}

Design SLHD::InnerGetOptimalSLHD(Design design, double p0, int max_iter_cnt) {
  auto opt_design = design.GetA();
  double opt_val = design.GetVal();
  double cur_val = opt_val;
  std::uniform_real_distribution<double> uniform_dis(std::nextafter(0.0, 1.0), 1.0);
  double temp = temp_init_;
  int cnt = 0;
  while (cnt < max_iter_cnt) {
    int i = 0;
    bool change_flag = false;
    while (i < imax_ && cnt < max_iter_cnt) {
      Log(cnt, design);
      ++cnt;

      int col = rng_() % k_;
      int r1, r2;
      double z = uniform_dis(rng_);
      if (z <= p0 || t_ == 1) {
        // Swap two elements in a slice
        int sid = rng_() % t_;
        r1 = rng_() % m_;
        r2 = Rand2(m_, r1);
        r1 = sid * m_ + r1;
        r2 = sid * m_ + r2;
      } else {
        // Swap two elements with same initial level between two slices
        int level = rng_() % m_;
        const auto& pos = design.GetPos()[col][level];
        r1 = rng_() % t_;
        r2 = Rand2(t_, r1);
        r1 = pos[r1];
        r2 = pos[r2];
      }

      design.SwapInCol(col, r1, r2);
      double tmp_val = design.GetVal();

      if (tmp_val < cur_val) {
        change_flag = true;
        cur_val = tmp_val;
      } else {
        double delta = tmp_val - cur_val;
        double acpt_prob = std::exp(-delta / temp);
        bool acpt = uniform_dis(rng_) <= acpt_prob;
        if (acpt) {
          change_flag = true;
          cur_val = tmp_val;
        } else {
          design.SwapInCol(col, r1, r2); // Restore
        }
      }

      if (cur_val < opt_val) {
        opt_val = cur_val;
        opt_design = design.GetA();
        i = 0;
      } else {
        ++i;
      }
    }
    if (!change_flag) break;
    temp *= (1 - temp_rate_);
  }
  return Design(opt_design, t_);
}
} // namespace SLHD