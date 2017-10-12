/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */




// calculation method for the cross sections
enum class ComputationMethod {Analytic, Lookup, Parametrized};

template <ComputationMethod method>
class PhotonCrossSection; 

template<>
class PhotonCrossSection <ComputationMethod::Analytic>
{
    using float_t = float;
public:
    float_t xs_pi_pi_rho0();
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

};

template<>
class PhotonCrossSection <ComputationMethod::Lookup>
{};

template<>
class PhotonCrossSection <ComputationMethod::Parametrized>
{};