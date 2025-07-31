#include <CLInterval.h>

#include <cmath>
#include <compare>
#include <format>
#include <iostream>

// Helper function for comparison operators
bool double_approx_equal(const double lhs, const double rhs) noexcept {
  static constexpr auto eps = 1e-20;    // For comparison with zero
  static constexpr auto relEps = 1e-5;  // For relative comparison between two numbers, one of which is nonzero
  return (std::isnan(lhs) && std::isnan(rhs)) || (std::abs(lhs) < eps && std::abs(rhs) < eps) ||
         std::abs((lhs - rhs) / std::max(std::abs(lhs), std::abs(rhs))) < relEps;
}

std::strong_ordering double_approx_three_way(const double lhs, const double rhs) noexcept {
  if (double_approx_equal(lhs, rhs)) return std::strong_ordering::equal;
  if (std::isnan(lhs) || lhs < rhs) return std::strong_ordering::less;
  return std::strong_ordering::greater;
}

bool CLInterval::operator==(const CLInterval& rhs) const noexcept {
  return this->minmethod == rhs.minmethod && this->maxmethod == rhs.maxmethod &&
         this->centralmethod == rhs.centralmethod && double_approx_equal(this->pvalue, rhs.pvalue) &&
         double_approx_equal(this->min, rhs.min) && double_approx_equal(this->max, rhs.max);
}

std::strong_ordering CLInterval::operator<=>(const CLInterval& rhs) const noexcept {
  if (auto cmp = this->minmethod <=> rhs.minmethod; cmp != 0) { return cmp; }
  if (auto cmp = this->maxmethod <=> rhs.maxmethod; cmp != 0) { return cmp; }
  if (auto cmp = this->centralmethod <=> rhs.centralmethod; cmp != 0) { return cmp; }
  if (auto cmp = double_approx_three_way(rhs.pvalue, this->pvalue); cmp != 0) { return cmp; }
  if (auto cmp = double_approx_three_way(this->min, rhs.min); cmp != 0) { return cmp; }
  if (auto cmp = double_approx_three_way(this->max, rhs.max); cmp != 0) { return cmp; }
  return std::strong_ordering::equal;
}

void CLInterval::print() const {
  std::cout << std::format("CLInterval {{\n"
                           "  p-value           = {:.4f}\n"
                           "  p-value at center = {:.4f}\n"
                           "  {:.4e} [{:.4e}, {:4e}]\n"
                           "  methods: {:s} [{:s}, {:s}]\n"
                           "  closed borders: [{:s}, {:s}]\n"
                           "}}",
                           pvalue, pvalueAtCentral, central, min, max, std::string(centralmethod),
                           std::string(minmethod), std::string(maxmethod), minclosed, maxclosed)
            << std::endl;
}
