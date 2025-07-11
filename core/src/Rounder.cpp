#include <Rounder.h>
#include <Utils.h>

using namespace Utils;

Rounder::Rounder(const OptParser* arg, double cllo, double clhi, double central)
    : m_cllo(cllo), m_clhi(clhi), m_central(central) {
  assert(arg);
  this->arg = arg;
}

int Rounder::getNsubdigits() const {
  if (arg->digits > -1) return arg->digits;
  return TMath::Max(calcNsubdigits(fabs(m_central - m_cllo)), calcNsubdigits(fabs(m_central - m_clhi)));
}

double Rounder::CLlo() const { return Round(m_cllo, getNsubdigits()); }

double Rounder::CLhi() const { return Round(m_clhi, getNsubdigits()); }

double Rounder::central() const { return Round(m_central, getNsubdigits()); }

///
/// Compute rounded negative error.
/// Use the rounded central values and interval boundaries for this.
/// Would we round the +/- errors themselves, we'd be geting inconsistent intervals!
/// \return minus error (always positive)
///
double Rounder::errNeg() const { return fabs(central() - CLlo()); }

double Rounder::errPos() const { return fabs(central() - CLhi()); }
