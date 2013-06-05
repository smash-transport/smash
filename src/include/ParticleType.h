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
 ParticleType() : mass_(-1), lifetime_(-1), isospin_(100), charge_(100),
    pdgcode_(-1), name_("unknown") {}
    /* Explicit constructor */
    ParticleType(std::string n, float m, int id) : mass_(m), lifetime_(-1),
    isospin_(1), charge_(0), pdgcode_(id), name_(n) {}
    /* set data */
    void inline set(const std::string &n, const float &m, const int &id,
                    const int &isospin, const int &charge);
    /* access data */
    int inline charge(void) const;
    int inline isospin(void) const;
    std::string inline name(void) const;
    float inline mass(void) const;
    int inline pdgcode(void) const;

  private:
    /* Data of the particle type */
    float mass_;
    float lifetime_;
    float isospin_;
    int charge_;
    int pdgcode_;
    std::string name_;
};

void inline ParticleType::set(const std::string &NAME, const float &MASS,
            const int &ID, const int &ISOSPIN, const int &CHARGE) {
  mass_ = MASS;
  pdgcode_ = ID;
  name_ = NAME;
  isospin_ = ISOSPIN;
  charge_ = CHARGE;
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

#endif  // SRC_INCLUDE_PARTICLETYPE_H_
