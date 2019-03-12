/*
 * File:   LocalExpansionCoeff.cpp
 * Author: matteo
 *
 * Created on April 3, 2015, 6:05 PM
 */

#include "local_expansion_coeff.h"
#include "string.h"
#include <iostream>

using std::cout;
using std::endl;

LocalExpansionCoeff::LocalExpansionCoeff()
{
  _p=0;
  _coeff= NULL;
}

LocalExpansionCoeff::LocalExpansionCoeff(const unsigned int &p)
{
  _p=p;
  _coeff= new double[(p+1)*(p+1)*(p+1)*(p+2)/2];
}

LocalExpansionCoeff::LocalExpansionCoeff(const LocalExpansionCoeff &orig)
{
  _p=orig._p;
  _coeff = new double[(_p+1)*(_p+1)*(_p+1)*(_p+2)/2];
  memcpy(_coeff, orig._coeff, this->getNumberOfElements());
}

LocalExpansionCoeff::~LocalExpansionCoeff()
{
  delete [] _coeff;
}

inline unsigned int LocalExpansionCoeff::getNumberOfElements()
{
  return (_p+1)*(_p+1)*(_p+1)*(_p+2)/2;
}

double LocalExpansionCoeff::get(const unsigned int &n, const unsigned int &m, const unsigned int &nn, const unsigned int &mm)
{

  return _coeff[ (mm+nn) +
                 getNNOffset(nn) +
                 getMOffset(m) +
                 getNOffset(n)];
}

void LocalExpansionCoeff::set(const unsigned int &n, const unsigned int &m, const unsigned int &nn, const unsigned int &mm, const double &value)
{

  _coeff[ (mm+nn) +
          getNNOffset(nn) +
          getMOffset(m) +
          getNOffset(n)]=value;
}

unsigned int LocalExpansionCoeff::getNNOffset(const unsigned int &nn)
{
  return (nn)*(nn);
}

unsigned int LocalExpansionCoeff::getMOffset(const unsigned int &m)
{
  return (_p+1)*(_p+1)*m;
}

unsigned int LocalExpansionCoeff::getNOffset(const unsigned int &n)
{
  return ((_p+1)*(_p+1)*(n+1)*(n))/2;
}


/*
 * Debugging and test of indexes
 */
//unsigned int const LocalExpansionCoeff::loopDebugger(const unsigned int & p){
//    unsigned int count=0;
//    for(unsigned int n=0; n< p+1; ++n){
//        for(unsigned int m=0; m < n+1; ++m){
//            for(int nn=0; nn< p+1; ++nn){
//                for(int mm=-nn; mm< nn+1; ++mm){
//                    ++count;
//                }
//            }
//        }
//    }
//    return count;
//}
//
//void LocalExpansionCoeff::fillCoeffWithIndex(){
//    unsigned int count=0;
//    for(unsigned int n=0; n< _p+1; ++n){
//        for(unsigned int m=0; m < n+1; ++m){
//            for(int nn=0; nn< _p+1; ++nn){
//                for(int mm=-nn; mm< nn+1; ++mm){
//                    ++count;
//                    _coeff[count-1]=count;
//
//                }
//            }
//        }
//    }
//}
//
//void LocalExpansionCoeff::printCoeff(){
//    unsigned int count=0;
//    for(unsigned int n=0; n< _p+1; ++n){
//        for(unsigned int m=0; m < n+1; ++m){
//            for(int nn=0; nn< _p+1; ++nn){
//                for(int mm=-nn; mm< nn+1; ++mm){
//                    ++count;
//                    cout << "(" << n << " " << m << " " << nn << " " << mm << ")"<< endl;
//                    cout << "_coeff: " << _coeff[count-1]<< endl;
//                    cout << "idx:" << (mm+nn) +
//                        getNNOffset(nn) +
//                        getMOffset(m)+
//                        getNOffset(n) << endl;
//                    cout << "get: " << this->get(n,m,nn,mm) << endl;
//                    cout << endl;
//
//                }
//            }
//        }
//    }
//}

