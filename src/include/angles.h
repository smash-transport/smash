/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ANGLES_H_
#define SRC_INCLUDE_ANGLES_H_

#include <iostream>
#include <stdexcept>

#include "constants.h"
#include "iomanipulators.h"
#include "random.h"
#include "threevector.h"

namespace smash {

/**
 * Angles provides a common interface for generating directions: i.e.,
 * two angles that should be interpreted as azimuthal and polar angles.
 *
 * Usage:
 * ------
 * \code
 * #include "Angles.h"
 *
 * Angles direction;
 *
 * direction.distribute_isotropically();
 * double azimuthal_angle = direction.phi();
 * double cosine_of_polar_angle = direction.costheta();
 * double x_projection_of_vector = direction.x();
 * direction.set_phi(0);
 * double new_azimuthal_angle = direction.phi();
 * // new_azimuthal_angle == 0.
 * \endcode
 *
 */
/*
 * Internals
 * ---------
 *
 * The object internally stores the azimuthal angle \f$\varphi\f$ and
 * the cosine of the polar angle \f$\cos\vartheta\f$. Nobody should rely
 * on this never changing, though; the interface user should be totally
 * oblivious to this.
 *
 * Possible future improvements
 * ----------------------------
 *
 * More distributions need to be implemented once there is a physics
 * case to use them.
 */
class Angles {
 public:
  /**
   * Default constructor.
   * \return Angles pointing in x-direction.
   */
  Angles() : phi_(0), costheta_(0) {}
  /**
   * \return Angles with given phi and cos(theta).
   * \param[in] ph the azimuthal angle
   * \param[in] cost cosine of the polar angle
   */
  Angles(double ph, double cost) {
    set_phi(ph);
    set_costheta(cost);
  }
  /**
   * Populate the object with a new direction.
   *
   * the direction is taken randomly from a homogeneous distribution,
   * i.e., each point on a unit sphere is equally likely.
   */
  void distribute_isotropically();
  /**
   * Sets the azimuthal angle.
   *
   * \param[in] phi Any real number to set the azimuthal angle \f$\varphi\f$
   * to.
   */
  void set_phi(const double phi);
  /**
   * Set the polar angle from its cosine.
   *
   * This is the preferred way of setting the polar information.
   *
   * \param[in] cos Cosine of the polar angle \f$\cos\vartheta\f$.
   * \throws InvalidTheta If cos is not in range [-1 .. 1].
   */
  void set_costheta(const double cos);
  /**
   * Set the polar angle.
   *
   * In the current implementation, costheta is stored inside the
   * object. Don't convert the angle to and from cosine, use the
   * set-function for the thing you have at hand.
   *
   * \param[in] theta Any real number to set the polar angle \f$\vartheta\f$
   * to.
   */
  void set_theta(const double theta);
  /**
   * Advance polar angle.
   *
   * A positive addition means that we go towards the southpole.
   *
   * \see add_to_theta(const double& delta, const bool& reverse)
   *
   * \param[in] delta Angle increment.
   * \return true if pole has been crossed.
   * \throws InvalidTheta If delta is not in [\f$-\pi\f$ .. \f$\pi\f$].
   */
  bool add_to_theta(const double delta);
  /**
   * Advance polar angle
   *
   * When crossing a pole, azimuthal angle is changed by 180 degrees.
   *
   * \param[in] delta Angle increment.
   * \param[in] reverse If true, we start in the "far" hemisphere, meaning a
   * positive delta will shift the object towards the north pole.
   *
   * \return true if we end up in the far hemisphere, false if
   * we end up in the original hemisphere.
   * \throws InvalidTheta If delta is not in [\f$-\pi\f$ .. \f$\pi\f$].

   */
  bool add_to_theta(const double delta, const bool reverse);
  /// \return Azimuthal angle.
  double phi() const;
  /// \return Cosine of polar angle.
  double costheta() const;
  /// \return Sine of polar angle.
  double sintheta() const;
  /**
   * \return \f$x\f$ projection of the direction.
   *
   * \f$x = \sin\vartheta \cos\varphi\f$
   */
  double x() const;
  /**
   * \return \f$y\f$ projection of the direction.
   *
   * \f$y = \sin\vartheta \sin\varphi\f$
   */
  double y() const;
  /**
   * \return \f$z\f$ projection of the direction.
   *
   * \f$z = \cos\vartheta\f$
   **/
  double z() const;
  /// \return The unit three-vector corresponding to the angles.
  ThreeVector inline threevec() const;
  /// \return The polar angle.
  double theta() const;

  /**
   * \ingroup exception
   * Thrown for invalid values for theta.
   */
  struct InvalidTheta : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };

 private:
  /// Azimuthal angle \f$\varphi\f$.
  double phi_;
  /// Cosine of polar angle \f$\cos\vartheta\f$.
  double costheta_;
};

/**
 * \ingroup logging
 * Creates output for an Angles object in the form "φ: 0.1294, cos ϑ:  0.423".
 */
inline std::ostream &operator<<(std::ostream &out, const Angles &a) {
  return out << "φ:" << field << a.phi() << ", cos ϑ:" << field << a.costheta();
}

void inline Angles::distribute_isotropically() {
  /* Isotropic distribution: phi in [0, 2pi) and cos(theta) in [-1,1]. */
  phi_ = random::uniform(0.0, twopi);
  costheta_ = random::uniform(-1.0, 1.0);
}

void inline Angles::set_phi(const double newphi) {
  /* Make sure that phi is in the range [0,2pi).  */
  phi_ = newphi;
  if (newphi < 0 || newphi >= twopi) {
    phi_ -= twopi * std::floor(newphi / twopi);
  }
}

void inline Angles::set_costheta(const double newcos) {
  costheta_ = newcos;
  /* check if costheta_ is in -1..1. If not, well. Error handling here
   * is a lot harder than in the above. Still, I will silently do the
   * same as above. Note, though, that costheta = 1 is allowed, even if
   * it cannot be generated by distribute_isotropically(). */
  if ((costheta_ < -1. - really_small) || (costheta_ > 1. + really_small)) {
    throw InvalidTheta("Wrong value for costheta (must be in [-1,1]): " +
                       std::to_string(costheta_));
  }
  if (costheta_ < -1.) {
    costheta_ = -1.;
  } else if (costheta_ > 1.) {
    costheta_ = 1.;
  }
}
void inline Angles::set_theta(const double newtheta) {
  /* no error handling necessary, because this gives a sensible answer
   * for every real number. */
  set_costheta(std::cos(newtheta));
}

bool inline Angles::add_to_theta(const double delta) {
  if (delta < -M_PI || delta > M_PI) {
    throw InvalidTheta("Cannot advance polar angle by " +
                       std::to_string(delta));
  }
  double theta_plus_delta = delta + theta();
  /* if sum is not in [0, PI], force it to be there:
   * "upper" overflow:
   * theta + delta + the_new_angle = 2*M_PI */
  if (theta_plus_delta > M_PI) {
    set_theta(twopi - theta_plus_delta);
    // set_phi takes care that phi_ is in [0 .. 2*M_PI]
    set_phi(phi() + M_PI);
    return true;  // meaning "we did change phi"
    // "lower" overflow: theta + delta switches sign
  } else if (theta_plus_delta < 0) {
    set_theta(-theta_plus_delta);
    set_phi(phi() + M_PI);
    return true;  // meaning "we did change phi"
    // no overflow: set theta, do not touch phi:
  } else {
    set_theta(theta_plus_delta);
  }
  return false;  // meaning "we did NOT change phi"
}
bool inline Angles::add_to_theta(const double delta, const bool reverse) {
  double plusminus_one = reverse ? -1.0 : +1.0;
  bool this_reverse = add_to_theta(plusminus_one * delta);
  /* if we had to reverse first time and now reverse again OR if we
   * didn't reverse in either part, we do not reverse in total.
   * else: if we reverse in one, but not the other part, we reverse in
   * total. */
  return this_reverse ^ reverse;
}

double inline Angles::costheta() const { return costheta_; }
double inline Angles::phi() const { return phi_; }
double inline Angles::sintheta() const {
  return std::sqrt(1.0 - costheta_ * costheta_);
}
double inline Angles::x() const { return sintheta() * std::cos(phi_); }
double inline Angles::y() const { return sintheta() * std::sin(phi_); }
double inline Angles::z() const { return costheta_; }

ThreeVector inline Angles::threevec() const {
  return ThreeVector(x(), y(), z());
}

double inline Angles::theta() const { return std::acos(costheta_); }

}  // namespace smash

#endif  // SRC_INCLUDE_ANGLES_H_
