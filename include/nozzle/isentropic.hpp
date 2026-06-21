// Copyright Roman Dial 2026

#pragma once

#include <cmath>
#include <cassert>
#include "numerics/root_finding.hpp"

namespace nozzle::core {

    /* Exit Mach from stagnation to exit pressure ratio
    
       p0_over_pe: p0 / pe  [-]
       gamma: ratio of specific heats [-]
       returns: exit Mach number [-]

    */

    inline double exit_mach_from_pressure_ratio(const double& p0_over_pe, const double& gamma) {
        const double exponent = (gamma - 1.0) / gamma;
        const double bracket  = std::pow(p0_over_pe, exponent) - 1.0;
        return std::sqrt((2.0 / (gamma - 1.0)) * bracket);
    }   

    /* Area Ratio Ae/A* from Mach number 
       
       returns: area ratio [-] 
   */

    inline double area_ratio_from_mach(const double& mach, const double& gamma) {
        const double gp1 = gamma + 1.0;
        const double gm1 = gamma - 1.0;
        const double inner = 1.0 + 0.5 * gm1 * mach * mach;
        const double exponent = gp1 / (2.0 * gm1);
        return (1.0 / mach) * std::pow((2.0 / gp1) * inner, exponent);
    }

    /* Prandtl-Meyer angle: the total turning angle a supersonic flow has
       undergone, expanding from Mach 1 to the given Mach number.
       
       returns:
        mach: local Mach number (>= 1) [-]
        gamma: ratio of specific heats [-]
        returns: Prandtl-Meyer angle nu [rad]
    */
    inline double prandtl_meyer(const double mach, const double gamma) {
        const double msm1 = mach * mach - 1.0;
        const double gp1 = gamma + 1.0;
        const double gm1 = gamma - 1.0;
        const double gratio = std::sqrt(gp1 / gm1);
        const double inner = std::sqrt((gm1 / gp1) * msm1);
        return gratio * std::atan(inner) - std::atan(std::sqrt(msm1));
    }

    /* Exit temperature from chamber temp and exit Mach
        
         T0: chamber temperature [K]
         Me: exit Mach number [-]
         gamma: ratio of specific heats [-]
         
         returns: exit temperature [K]
    */

    inline double exit_temperature(const double& T0, const double& mach, const double& gamma) {
        const double inner = 1.0 + 0.5 * (gamma - 1.0) * mach * mach;
        return T0 / inner;
    }

    /* Exit Velocity from chamber temp and exit Mach
        
         T0: chamber temperature [K]
         Me: exit Mach number [-]
         gamma: ratio of specific heats [-]
         R: specific gas constant [J/(kg·K)]
         
         returns: exit velocity [m/s]
    */

    inline double exit_velocity(const double& T0, const double& mach, const double & gamma, const double& R) {
        const double Te = exit_temperature(T0, mach, gamma);
        const double sound_speed = std::sqrt(gamma * R * Te);
        return mach * sound_speed;
    }

    /* Specific impulse */

    inline double specific_impulse(const double& exit_velocity) {
        constexpr double g0 = 9.80665; // standard gravity [m/s^2]
        return exit_velocity / g0;
    }

    /* Thrust Coefficient */

    inline double thrust_coefficient(const double gamma, const double pe_over_p0, const double pa_over_p0, const double area_ratio) {
        const double gp1 = gamma + 1.0;
        const double gm1 = gamma - 1.0;
        const double gterm1 = (2.0 * gamma * gamma) / gm1;
        const double gterm2 = 2.0 / gp1; 
        const double pimterm = std::pow(pe_over_p0, (gm1 / gamma));
        const double momentum_term = std::sqrt(gterm1 * std::pow(gterm2, (gp1 / gm1)) * (1.0 - pimterm));
        const double pressure_term = (pe_over_p0 - pa_over_p0) * area_ratio;
        return momentum_term + pressure_term;
    }

    /* Mach from Prandtl Meyer 
    
       returns: Mach [-] 
    */
    inline double mach_from_prandtl_meyer(double nu_target, double gamma) {
        
        auto f = [gamma, nu_target](double M) {
            return prandtl_meyer(M, gamma) - nu_target;
        };
        
        auto df = [gamma](double M) {
            const double msm1 = M * M - 1.0;
            const double isofactor = 1.0 + 0.5 * (gamma - 1.0) * M * M;
            return std::sqrt(msm1) / (M * isofactor);
        };

        return numerics::safeguarded_newton(f, df, 1.0, 100.0);
    }

    /* Mach angle
    
       returns: Mach angle [rad]
    */
    inline double mach_angle(double mach) {
        assert(mach >= 1.0);
        return std::asin(1.0 / mach);
    }

}   // namespace nozzle::core

