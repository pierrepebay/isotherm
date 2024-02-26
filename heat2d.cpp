#include <iostream>
#include <array>

const int N = 5;

uint64_t cartesian_to_index(uint64_t i, uint64_t j) {
  return j * N + i;
}

int main() {
  std::array<double, N*N> T{};

  for (double & u : T)
    u = 0.0;

  for (uint64_t i = 0; i < N; ++i) {
    for (uint64_t j = 0; j < N; j++) {
      auto index = cartesian_to_index(i, j);
      std::cout << "Cartesian: (" << i << ", " << j << ")" << " -- Index: " << index << std::endl;
    }
  }
}
