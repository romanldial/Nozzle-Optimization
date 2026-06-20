// Copyright Roman Dial 2026

#pragma once

#include <cmath>

namespace nozzle::core {

    /* Exit Mach from stagnation to exit pressure ratio
    
       p0_over_pe: p0 / pe  [-]
       gamma: ratio of specific heats [-]
       returns: exit Mach number [-]

    */

    inline double exit_mach_from_pressure_ratio(const double& p0_over_pe,
                                                const double& gamma) {
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


}   // namespace nozzle::core

