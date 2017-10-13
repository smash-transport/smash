/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <cmath>


// calculation method for the cross sections
enum class ComputationMethod {Analytic, Lookup, Parametrized};


// usage would be 
// PhotonCrossSection xs<Analytic>; xs.xs_pi_pi_rho(var1, var2...)
template <ComputationMethod method>
class PhotonCrossSection; 

template<>
class PhotonCrossSection <ComputationMethod::Analytic>
{
    
public: 
    using float_t = float;

    // could be static 
    // TODO: Clean up parameters. Some are unused. 
    float_t xs_pi_pi_rho0(const float_t m1, const float_t m2, 
        const float_t m3, const float_t t1, const float_t t2, const float_t s,
        const float_t mpion, const float_t mrho);
    float_t xs_pi0_pi_rho(const float_t m1, const float_t m2, 
        const float_t m3, const float_t t1, const float_t t2, const float_t s,
        const float_t mpion, const float_t mrho);
    float_t xs_pi_rho0_pi(const float_t m1, const float_t m2, 
        const float_t m3, const float_t t1, const float_t t2, const float_t s,
        const float_t mpion, const float_t mrho);
    float_t xs_pi_rho_pi0(const float_t m1, const float_t m2, 
        const float_t m3, const float_t t1, const float_t t2, const float_t s,
        const float_t mpion, const float_t mrho);
    float_t xs_pi0_rho_pi(const float_t m1, const float_t m2, 
        const float_t m3, const float_t t1, const float_t t2, const float_t s,
        const float_t mpion, const float_t mrho);
    float_t xs_pi0_rho0_pi0(const float_t m1, const float_t m2, 
        const float_t m3, const float_t t1, const float_t t2, const float_t s,
        const float_t mpion, const float_t mrho);
    
    
        // and so on
    
private:

    const float_t to_mb = 0.3894;
    const float_t Const = 0.059;
    const float_t g_POR = 11.93;
    const float_t ma1 = 1.26;
    const float_t ghat = 6.4483;
    const float_t eta1 = 2.3920;
    const float_t eta2 = 1.9430;
    const float_t delta = -0.6426;
    const float_t C4 = -0.14095;
    const float_t Gammaa1 = 0.4;
    const float_t Pi = M_PI;
    const float_t m_omega = 0.783;
    

};

template<>
class PhotonCrossSection <ComputationMethod::Lookup>
{};

template<>
class PhotonCrossSection <ComputationMethod::Parametrized>
{};