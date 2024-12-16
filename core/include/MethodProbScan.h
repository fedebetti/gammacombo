/*
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: August 2012
 *
 */

#ifndef MethodProb_h
#define MethodProb_h

#include <TH1F.h>

#include "MethodAbsScan.h"

class MethodProbScan : public MethodAbsScan
{
public:
    MethodProbScan(Combiner* comb);
    MethodProbScan(OptParser* opt);
    MethodProbScan();
    virtual ~MethodProbScan();

    virtual int   computeCLvalues();  // compute CL histograms depending on desired test statistic
    float         getChi2min(float scanpoint) const;
    inline TH1F*  getHChi2min() const {return hChi2min;};
    void          saveSolutions();
    void          saveSolutions2d();
    virtual int   scan1d(bool fast=false, bool reverse=false, bool quiet=false);
    virtual int   scan2d();
    inline void   setScanDisableDragMode(bool f=true){scanDisableDragMode = f;};

protected:
    bool  computeInnerTurnCoords(const int iStart, const int jStart, const int i, const int j,
                                 int &iResult, int &jResult, int nTurn) const;
    bool  deleteIfNotInCurveResults2d(RooSlimFitResult *r);
    void  sanityChecks() const;
    bool  scanDisableDragMode = false;
    int   nScansDone = 0;  // count the number of times a scan was done
};

#endif
