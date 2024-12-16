/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: April 2013
 *
 * Class to round confidence intervals and central values.
 *
 **/

#ifndef Rounder_h
#define Rounder_h

#include "OptParser.h"
#include "Utils.h"

class Rounder
{
public:

    Rounder(const OptParser *arg, float cllo, float clhi, float central);
    ~Rounder();
    
    int   getNsubdigits() const;
    float CLlo() const;
    float CLhi() const;
    float central() const;
    float errNeg() const;
    float errPos() const;
    
private:
      
    const OptParser *arg;  ///< command line arguments
    float m_cllo;
    float m_clhi;
    float m_central;
};

#endif
