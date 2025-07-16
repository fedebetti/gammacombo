#include <CLInterval.h>

#include <format>
#include <iostream>

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
