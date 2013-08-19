/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
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
    void inline set(const std::string &n, const float &m, const float &w,
                    const int &id, const int &isosp, const int &ch,
                    const int &sp);
    /* access data */
    std::string inline name(void) const;
    float inline mass(void) const;
    float inline width(void) const;
    int inline pdgcode(void) const;
    /* Isospin is 2 * particle data book value */
    int inline isospin(void) const;
    int inline charge(void) const;
    /* Spin is 2 * particle data book value */
    int inline spin(void) const;

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

void inline ParticleType::set(const std::string &NAME, const float &MASS,
     const float &WIDTH, const int &ID, const int &ISOSPIN, const int &CHARGE,
     const int &SPIN) {
  mass_ = MASS;
  width_ = WIDTH;
  pdgcode_ = ID;
  name_ = NAME;
  isospin_ = ISOSPIN;
  charge_ = CHARGE;
  spin_ = SPIN;
}

int inline ParticleType::charge(void) const {
  return charge_;
}

int inline ParticleType::isospin(void) const {
  return isospin_;
}

float inline ParticleType::mass(void) const {
  return mass_;
}

std::string inline ParticleType::name(void) const {
  return name_;
}

int inline ParticleType::pdgcode(void) const {
  return pdgcode_;
}

int inline ParticleType::spin(void) const {
  return spin_;
}

float inline ParticleType::width(void) const {
  return width_;
}

#endif  // SRC_INCLUDE_PARTICLETYPE_H_
