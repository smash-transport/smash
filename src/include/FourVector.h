/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_FOURVECTOR_H_
#define SRC_INCLUDE_FOURVECTOR_H_

#include <cmath>

class FourVector {
 public:
    /* default constructor */
    FourVector(): x0_(0.0), x1_(0.0), x2_(0.0), x3_(0.0) {}
    /* useful constructor */
    FourVector(double y0, double y1, double y2, double y3): x0_(y0),
      x1_(y1), x2_(y2), x3_(y3) {}
    /* t, z, x_\perp */
    double inline x0(void) const;
    void inline set_x0(double t);
    double inline x1(void) const;
    void inline set_x1(double z);
    double inline x2(void) const;
    void inline set_x2(double x);
    double inline x3(void) const;
    void inline set_x3(double y);
    /* set all four values */
    void inline set_FourVector(const double t, const double z, const double x,
      const double y);
    /* inlined operations */
    double inline Dot(const FourVector &a) const;
    double inline Dot() const;
    double inline DotThree(const FourVector &a) const;
    double inline DotThree() const;
    double inline DiffThree(const FourVector &a) const;
    /* operations */
    FourVector LorentzBoost(const FourVector &b) const;

    /* overloaded operators */
    bool inline operator==(const FourVector &a) const;
    bool inline operator!=(const FourVector &a) const;
    bool inline operator<(const FourVector &a) const;
    bool inline operator>(const FourVector &a) const;
    bool inline operator<=(const FourVector &a) const;
    bool inline operator>=(const FourVector &a) const;
    bool inline operator==(const double &a) const;
    bool inline operator!=(const double &a) const;
    bool inline operator<(const double &a) const;
    bool inline operator>(const double &a) const;
    bool inline operator<=(const double &a) const;
    bool inline operator>=(const double &a) const;
    FourVector inline operator+=(const FourVector &a);
    FourVector inline operator-=(const FourVector &a);
    FourVector inline operator*=(const double &a);
    FourVector inline operator/=(const double &a);

 private:
    double x0_, x1_, x2_, x3_;
};

double inline FourVector::x0(void) const {
  return x0_;
}

void inline FourVector::set_x0(const double t) {
  x0_ = t;
}

double inline FourVector::x1(void) const {
  return x1_;
}

void inline FourVector::set_x1(const double z) {
  x1_ = z;
}

double inline FourVector::x2(void) const {
  return x2_;
}

void inline FourVector::set_x2(const double x) {
  x2_ = x;
}

double inline FourVector::x3(void) const {
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

/* all four vector components are equal */
bool inline FourVector::operator==(const FourVector &a) const {
  return fabs(x0_ - a.x0_) < 1e-12 && fabs(x1_ - a.x1_) < 1e-12
    && fabs(x2_ - a.x2_) < 1e-12 && fabs(x3_ - a.x3_) < 1e-12;
}

/* use == operator for the inverse */
bool inline FourVector::operator!=(const FourVector &a) const {
  return !(*this == a);
}

/* all four vector components are below comparison vector */
bool inline FourVector::operator<(const FourVector &a) const {
  return (x0_ < a.x0_) && (x1_ < a.x1_) && (x2_ < a.x2_) && (x3_ < a.x3_);
}

/* use < operator for the inverse by switching arguments */
bool inline FourVector::operator>(const FourVector &a) const {
  return a < *this;
}

/* use > operator for less equal */
bool inline FourVector::operator<=(const FourVector &a) const {
  return !(*this > a);
}

/* use < operator for greater equal */
bool inline FourVector::operator>=(const FourVector &a) const {
  return !(*this < a);
}

/* all vector components are equal to that number */
bool inline FourVector::operator==(const double &a) const {
  return fabs(x0_ - a) < 1e-12 && fabs(x1_ - a) < 1e-12
    && fabs(x2_ - a) < 1e-12 && fabs(x3_ - a) < 1e-12;
}

/* use == operator for the inverse */
bool inline FourVector::operator!=(const double &a) const {
  return !(*this == a);
}

/* all vector components are below that number */
bool inline FourVector::operator<(const double &a) const {
  return (x0_ < a) && (x1_ < a) && (x2_ < a) && (x3_ < a);
}

/* all vector components are above that number */
bool inline FourVector::operator>(const double &a) const {
  return (x0_ > a) && (x1_ > a) && (x2_ > a) && (x3_ > a);
}

/* all vector components are less equal that number */
bool inline FourVector::operator<=(const double &a) const {
  return !(*this > a);
}

/* all vector components are greater equal that number */
bool inline FourVector::operator>=(const double &a) const {
  return !(*this < a);
}

/* assignement addition */
FourVector inline FourVector::operator+=(const FourVector &a) {
  this->x0_ += a.x0_;
  this->x1_ += a.x1_;
  this->x2_ += a.x2_;
  this->x3_ += a.x3_;
  return *this;
}

/* addition uses += */
inline FourVector operator+(FourVector a, const FourVector &b) {
  a += b;
  return a;
}

/* assignement subtraction */
FourVector inline FourVector::operator-=(const FourVector &a) {
  this->x0_ -= a.x0_;
  this->x1_ -= a.x1_;
  this->x2_ -= a.x2_;
  this->x3_ -= a.x3_;
  return *this;
}

/* subtraction uses -= */
inline FourVector operator-(FourVector a, const FourVector &b) {
  a -= b;
  return a;
}

/* assignement factor multiplication */
FourVector inline FourVector::operator*=(const double &a) {
  this->x0_ *= a;
  this->x1_ *= a;
  this->x2_ *= a;
  this->x3_ *= a;
  return *this;
}

/* factor multiplication uses *= */
inline FourVector operator*(FourVector a, const double &b) {
  a *= b;
  return a;
}

/* assignement factor division */
FourVector inline FourVector::operator/=(const double &a) {
  this->x0_ /= a;
  this->x1_ /= a;
  this->x2_ /= a;
  this->x3_ /= a;
  return *this;
}

/* factor division uses /= */
inline FourVector operator/(FourVector a, const double &b) {
  a /= b;
  return a;
}

double inline FourVector::Dot(const FourVector &a) const {
  return x0_ * a.x0_ - x1_ * a.x1_ - x2_ * a.x2_ - x3_ * a.x3_;
}

double inline FourVector::Dot() const {
  return x0_ * x0_ - x1_ * x1_ - x2_ * x2_ - x3_ * x3_;
}

double inline FourVector::DotThree(const FourVector &a) const {
  return - x1_ * a.x1_ - x2_ * a.x2_ - x3_ * a.x3_;
}

double inline FourVector::DotThree() const {
  return - x1_ * x1_ - x2_ * x2_ - x3_ * x3_;
}

double inline FourVector::DiffThree(const FourVector &a) const {
  return x1_ - a.x1_ + x2_ - a.x2_ + x3_ - a.x3_;
}

#endif  // SRC_INCLUDE_FOURVECTOR_H_
