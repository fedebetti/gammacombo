#include <Rounder.h>
#include <Utils.h>

#include <cmath>

using Utils::Round;

Rounder::Rounder(const OptParser* arg, const double clMin, double clMax, double central)
    : m_clMin(clMin), m_clMax(clMax), m_central(central) {
  assert(arg);
  this->arg = arg;
  if (this->arg->digits > -1) {
    m_nSubDigits = arg->digits;
  } else {
    // Base the number of digits to show on the parameters uncertainty or, if the central value is not available,
    // on the width of the CL interval.
    using Utils::calcNsubdigits;
    if (!std::isnan(m_central))
      m_nSubDigits =
          std::max(calcNsubdigits(std::abs(m_central - m_clMin)), calcNsubdigits(std::abs(m_central - m_clMax)));
    else
      m_nSubDigits = calcNsubdigits(m_clMax - m_clMin);
  }
}

double Rounder::central() const { return Round(m_central, m_nSubDigits); }

double Rounder::errNeg() const { return Round(std::abs(m_central - m_clMin), m_nSubDigits); }

double Rounder::errPos() const { return Round(std::abs(m_central - m_clMax), m_nSubDigits); }

double Rounder::CLlo() const {
  const auto clMin = central() - errNeg();
  return (std::isnan(clMin)) ? Round(m_clMin, m_nSubDigits) : clMin;
}

double Rounder::CLhi() const {
  const auto clMax = central() + errPos();
  return (std::isnan(clMax)) ? Round(m_clMax, m_nSubDigits) : clMax;
}
