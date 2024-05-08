#include <vector>
#include <iostream>
#include <array>

namespace {
  struct dictionary {
    std::vector<std::array<double, 10> > vxyzp;
    std::vector<std::vector<std::vector<double>>> vx;
  };
}
