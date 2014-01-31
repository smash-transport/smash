/*
 *    Copyright (c) 2012-2013
 *      SMASH Team
 * 
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PARTICLETYPE_H_
#define SRC_INCLUDE_PARTICLETYPE_H_

#include <string>

class ParticleType {
 public:
    /* Use improbable values for default constructor */
    ParticleType() : name_("unknown"), mass_(-1), width_(-1), pdgcode_(-1),
                     isospin_(100), charge_(100), spin_(100) {}
    /* Explicit constructor */
    ParticleType(std::string n, float m, float w, int id, int isosp,
                 int ch, int sp) : name_(n), mass_(m), width_(w), pdgcode_(id),
                 isospin_(isosp), charge_(ch), spin_(sp) {}
    /* set data */
    inline void set(const std::string &n, float m, float w,
                    int id, int isosp, int ch, int sp);
    /* access data */
    inline std::string name(void) const;
    inline float mass(void) const;
    inline float width(void) const;
    inline int pdgcode(void) const;
    /* Isospin is 2 * particle data book value */
    inline int isospin(void) const;
    inline int charge(void) const;
    /* Spin is 2 * particle data book value */
    inline int spin(void) const;

 private:
    /* Data of the particle type */
    std::string name_;
    float mass_;
    float width_;
    int pdgcode_;
    int isospin_;
    int charge_;
    int spin_;
};

inline void ParticleType::set(const std::string &NAME, float MASS,
     float WIDTH, int ID, int ISOSPIN, int CHARGE, int SPIN) {
  mass_ = MASS;
  width_ = WIDTH;
  pdgcode_ = ID;
  name_ = NAME;
  isospin_ = ISOSPIN;
  charge_ = CHARGE;
  spin_ = SPIN;
}

inline int ParticleType::charge(void) const {
  return charge_;
}

inline int ParticleType::isospin(void) const {
  return isospin_;
}

inline float ParticleType::mass(void) const {
  return mass_;
}

inline std::string ParticleType::name(void) const {
  return name_;
}

inline int ParticleType::pdgcode(void) const {
  return pdgcode_;
}

inline int ParticleType::spin(void) const {
  return spin_;
}

inline float ParticleType::width(void) const {
  return width_;
}

#endif  // SRC_INCLUDE_PARTICLETYPE_H_
