/*
 *
 *    Copyright (c) 2014
 *      Bjørn Bäuchle <baeuchle@fias.uni-frankfurt.de>
 *      
 *    GNU General Public License (GPLv3)
 *
 */
 
 
#ifndef ANGLES_H
#define ANGLES_H

#include <cmath>


class angles {
 public:
   angles();
   // get new angles:
   void distribute_isotropously();
   // update angles:
   void set_phi(const double& phi);
   // set cos(theta) (preferred):
   void set_cos(const double& cos);
   // set theta (if you only have the angle, not the cosine):
   void set_theta(const double& theta);
   // get elements:
   double phi() const;
   double costheta() const;
   // sqrt(1-costheta**2)
   double sintheta() const;
   // x, y and z together give a normalized three vector which, thus,
   // has to be multiplied by its length.
   double x() const;
   double y() const;
   double z() const;
   // in case we really need the polar angle instead of just the cosine,
   // we can return acos(costheta).
   double theta() const;

 private:
   double phi_;
   double costheta_;
};

inline angles::angles() {}

void inline angles::distribute_isotropously() {
  // isotropous distribution: phi in [0, 2pi) and cos(theta) in [-1,1]
  phi_ = 2.0 * M_PI * drand48();
  costheta_ = -1.0 + 2.0 * drand48();
}

void inline angles::set_phi (const double& newphi) {
  phi_ = newphi;
}
void inline angles::set_cos (const double& newcos) {
  costheta_ = newcos;
}
void inline angles::set_theta (const double& newtheta) {
  set_cos( cos(newtheta) );
}

double inline angles::costheta() const { return costheta_; }
double inline angles::phi() const { return phi_; }
double inline angles::sintheta() const { return sqrt(1.0 - costheta_*costheta_); }
double inline angles::x() const { return costheta_*cos(phi_); }
double inline angles::y() const { return costheta_*sin(phi_); }
double inline angles::z() const { return sintheta(); }
double inline angles::theta() const { return acos(costheta_); }

#endif // ANGLES_H
