#ifndef SRC_INCLUDE_SMASH_PARTICLEQUANTITY_H_
#define SRC_INCLUDE_SMASH_PARTICLEQUANTITY_H_

#include "particledata.h"
#include <variant>
#include <functional>
#include <string>
#include <map>

namespace smash {

    
 class ParticleQuantity {
    public:
        using returnVariant = std::variant<int, double>;
        using FuncType = std::function<returnVariant(const particledata&)>;

        
        template<typename T>
        static FuncType createFunc(T (particledata::*memberfunc)() const) {
            return [memberfunc](const particledata& data) -> returnVariant {
                return (data.*memberfunc)();
            };
        }
        

        const std::string getUnit() const { return _unit; }
        const std::string getQuantity() const { return _quantity; }

        FuncType getValue;

        const static std::map<std::string, ParticleQuantity> ParticleQuantities;

    private:
        ParticleQuantity(std::string quantity, std::string unit, FuncType getValueFunc)
            : _quantity(std::move(quantity)), _unit(std::move(unit)), getValue(std::move(getValueFunc)) {}

        std::string _unit;
        std::string _quantity;
    };

    // Static member initialization outside the class definition
    const std::map<std::string, ParticleQuantity> ParticleQuantity::ParticleQuantities = {
        {"formation_time", {"formation_time", "fm", ParticleQuantity::createFunc(&particledata::formation_time)}}
    };

}

#endif  // SRC_INCLUDE_SMASH_PARTICLEQUANTITY_H_
