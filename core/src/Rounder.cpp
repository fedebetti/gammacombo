#include <Rounder.h>
#include <Utils.h>

using namespace Utils;

Rounder::Rounder(const OptParser* arg, float cllo, float clhi, float central)
    : m_cllo(cllo), m_clhi(clhi), m_central(central) {
  assert(arg);
  this->arg = arg;
}

Rounder::~Rounder() {}

int Rounder::getNsubdigits() const {
  if (arg->digits > -1) return arg->digits;
  return TMath::Max(calcNsubdigits(fabs(m_central - m_cllo)), calcNsubdigits(fabs(m_central - m_clhi)));
}

float Rounder::CLlo() const { return Round(m_cllo, getNsubdigits()); }

float Rounder::CLhi() const { return Round(m_clhi, getNsubdigits()); }

float Rounder::central() const { return Round(m_central, getNsubdigits()); }

///
/// Compute rounded negative error.
/// Use the rounded central values and interval boundaries for this.
/// Would we round the +/- errors themselves, we'd be geting inconsistent intervals!
/// \return minus error (always positive)
///
float Rounder::errNeg() const { return fabs(central() - CLlo()); }

float Rounder::errPos() const { return fabs(central() - CLhi()); }
