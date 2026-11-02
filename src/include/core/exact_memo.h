// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// exact_memo.h — the last results of a pure function, keyed by the bits of its inputs
//
// ExactMemo<NK, NV, K> keeps the K most recent (input, output) pairs of a function
// whose NV outputs are determined by its NK double inputs alone.  A lookup compares
// the inputs bit for bit, so a hit returns exactly the outputs the function returns
// for those inputs: a call answered from the memo and a computed call are identical,
// whatever the order of the calls.  Callers keep one instance per thread
// (thread_local), so no state is shared between threads.
// ---------------------------------------------------------------------------
#include <array>
#include <cstring>

namespace arche {

template <int NK, int NV, int K>
class ExactMemo {
 public:
  bool find(const std::array<double, NK>& in,
            std::array<double, NV>& out) const {
    for (int i = 0; i < n_; ++i) {
      if (std::memcmp(in_[i].data(), in.data(), sizeof(double) * NK) == 0) {
        out = out_[i];
        return true;
      }
    }
    return false;
  }

  // Replaces the oldest entry once all K are in use.
  void insert(const std::array<double, NK>& in,
              const std::array<double, NV>& out) {
    in_[next_] = in;
    out_[next_] = out;
    next_ = (next_ + 1) % K;
    if (n_ < K) ++n_;
  }

 private:
  std::array<std::array<double, NK>, K> in_{};
  std::array<std::array<double, NV>, K> out_{};
  int n_ = 0;
  int next_ = 0;
};

}  // namespace arche
