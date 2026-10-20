/*
 * File:   LocalExpansionCoeff.cpp
 * Author: matteo
 *
 * Created on April 3, 2015, 6:05 PM
 */

#include "local_expansion_coeff.h"

#include <iostream>

#include "string.h"

using std::cout;
using std::endl;

LocalExpansionCoeff::LocalExpansionCoeff()
{
  _p     = 0;
  _coeff = NULL;
}

LocalExpansionCoeff::LocalExpansionCoeff(const unsigned int &p)
{
  _p     = p;
  _coeff = new double[getNumberOfElements()];
}

LocalExpansionCoeff::LocalExpansionCoeff(const LocalExpansionCoeff &orig)
{
  _p     = orig._p;
  _coeff = new double[getNumberOfElements()];
  memcpy(_coeff, orig._coeff, this->getNumberOfElements());
}

LocalExpansionCoeff::~LocalExpansionCoeff()
{
  delete[] _coeff;
}

inline unsigned int
LocalExpansionCoeff::getNumberOfElements()
{
  return (_p + 1) * (_p + 1) * (_p + 1) * (_p + 2) / 2;
}

double
LocalExpansionCoeff::get(const unsigned int &n,
                         const unsigned int &m,
                         const unsigned int &nn,
                         const unsigned int &mm)
{
  return _coeff[(mm + nn) + getNNOffset(nn) + getMOffset(m) + getNOffset(n)];
}

void
LocalExpansionCoeff::set(const unsigned int &n,
                         const unsigned int &m,
                         const unsigned int &nn,
                         const unsigned int &mm,
                         const double       &value)
{
  _coeff[(mm + nn) + getNNOffset(nn) + getMOffset(m) + getNOffset(n)] = value;
}

unsigned int
LocalExpansionCoeff::getNNOffset(const unsigned int &nn)
{
  return (nn) * (nn);
}

unsigned int
LocalExpansionCoeff::getMOffset(const unsigned int &m)
{
  return (_p + 1) * (_p + 1) * m;
}

unsigned int
LocalExpansionCoeff::getNOffset(const unsigned int &n)
{
  return ((_p + 1) * (_p + 1) * (n + 1) * (n)) / 2;
}
