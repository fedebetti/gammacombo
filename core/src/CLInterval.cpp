#include "CLInterval.h"

using namespace std;

void CLInterval::print() const
{
    cout << "pvalue=" << pvalue
        << " pvalueAtCentral=" << pvalueAtCentral
        << " min=" << min
        << " max=" << max
        << " central=" << central
        << " minclosed=" << minclosed
        << " maxclosed=" << maxclosed
        << " minmethod=" << minmethod
        << " maxmethod=" << maxmethod
        << " centralmethod=" << centralmethod
        << endl;
}
