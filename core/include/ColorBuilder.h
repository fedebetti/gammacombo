/**
 * \author Till Moritz Karbach, moritz.karbach@cern.ch
 * \date August 2014
 * \brief Class that can make new darker or lighter colors.
 *
 **/

#ifndef ColorBuilder_h
#define ColorBuilder_h

class ColorBuilder {
 public:
  int darkcolor(int n) const;
  int darklightcolor(int n, double scale) const;
  int lightcolor(int n) const;
};

#endif
