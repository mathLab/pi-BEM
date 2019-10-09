// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2019 by the pi-BEM authors.
//
// This file is part of the pi-BEM library.
//
// The BEMStokes is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License version 2.1 as published by the Free Software Foundation.
// The full text of the license can be found in the file LICENSE at
// the top level of the BEMStokes distribution.
//
// Authors: Nicola Giuliani, Andrea Mola, Luca Heltai
//
// ---------------------------------------------------------------------

#ifndef LOCALEXPANSIONCOEFF_H
#define LOCALEXPANSIONCOEFF_H

class LocalExpansionCoeff
{
public:
  LocalExpansionCoeff();
  LocalExpansionCoeff(const unsigned int &p);
  LocalExpansionCoeff(const LocalExpansionCoeff &orig);
  double
  get(const unsigned int &n,
      const unsigned int &m,
      const unsigned int &nn,
      const unsigned int &mm);
  void
  set(const unsigned int &n,
      const unsigned int &m,
      const unsigned int &nn,
      const unsigned int &mm,
      const double &      value);
  unsigned int
  getNumberOfElements();
  unsigned int
  getNNOffset(const unsigned int &nn);
  unsigned int
  getMOffset(const unsigned int &m);
  unsigned int
  getNOffset(const unsigned int &n);
  virtual ~LocalExpansionCoeff();

  // Debugging and test of indexes
  //    static unsigned int const loopDebugger(const unsigned int & p);
  //    void fillCoeffWithIndex();
  //    void printCoeff();

private:
  unsigned int _p;
  double *     _coeff;
};

#endif /* LOCALEXPANSIONCOEFF_HPP */
