/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_FOURVECTOR_H_
#define SRC_INCLUDE_FOURVECTOR_H_

class FourVector {
  public:
    /* default constructor */
    FourVector(): x0_(0), x1_(0), x2_(0), x3_(0) {}
    FourVector(double y0, double y1, double y2, double y3): x0_(y0),
      x1_(y1), x2_(y2), x3_(y3) {}
    /* t, z, x_\perp */
    double inline x0(void);
    void inline set_x0(double t);
    double inline x1(void);
    void inline set_x1(double z);
    double inline x2(void);
    void inline set_x2(double x);
    double inline x3(void);
    void inline set_x3(double y);
    void inline set_FourVector(const double t, const double z, const double x,
      const double y);
    double Dot(FourVector);

    /* overloaded operators */
    FourVector inline operator+=(const FourVector &a);
    FourVector inline operator-=(const FourVector &a);
    FourVector inline operator*=(const double &a);

  private:
    double x0_, x1_, x2_, x3_;
};

double inline FourVector::x0(void) {
  return x0_;
}

void inline FourVector::set_x0(const double t) {
  x0_ = t;
}

double inline FourVector::x1(void) {
  return x1_;
}

void inline FourVector::set_x1(const double z) {
  x1_ = z;
}

double inline FourVector::x2(void) {
  return x2_;
}

void inline FourVector::set_x2(const double x) {
  x2_ = x;
}

double inline FourVector::x3(void) {
  return x3_;
}

void inline FourVector::set_x3(const double y) {
  x3_ = y;
}

void inline FourVector::set_FourVector(const double t, const double z,
                                       const double x, const double y) {
  x0_ = t;
  x1_ = z;
  x2_ = x;
  x3_ = y;
}

FourVector inline FourVector::operator+=(const FourVector &a) {
  this->x0_ += a.x0_;
  this->x1_ += a.x1_;
  this->x2_ += a.x2_;
  this->x3_ += a.x3_;
  return *this;
}

inline FourVector operator+(FourVector a, const FourVector &b) {
  a += b;
  return a;
}

FourVector inline FourVector::operator-=(const FourVector &a) {
  this->x0_ -= a.x0_;
  this->x1_ -= a.x1_;
  this->x2_ -= a.x2_;
  this->x3_ -= a.x3_;
  return *this;
}

inline FourVector operator-(FourVector a, const FourVector &b) {
  a -= b;
  return a;
}

FourVector inline FourVector::operator*=(const double &a) {
  this->x0_ *= a;
  this->x1_ *= a;
  this->x2_ *= a;
  this->x3_ *= a;
  return *this;
}

#endif  // SRC_INCLUDE_FOURVECTOR_H_
