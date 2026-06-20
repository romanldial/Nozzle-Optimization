// Copyright Roman Dial 2026

#pragma once

#include "chamber_conditions.hpp"
#include "design_result.hpp"
#include "isentropic.hpp"

namespace nozzle::design {

    /* Design Modes for the Nozzle Optimization 
        
        Target: pe = pambient 
        
    */

    inline DesignResult design_optimal_expansion(const ChamberConditions& chamber, const double& p_ambient) {
        const double p0_over_pe = chamber.p0 / p_ambient;
        const double Me = core::exit_mach_from_pressure_ratio(p0_over_pe, chamber.gamma);
        const double ar = core::area_ratio_from_mach(Me, chamber.gamma);
        const double ve = core::exit_velocity(chamber.T0, Me, chamber.gamma, chamber.R);
        const double isp = core::specific_impulse(ve);
        
        DesignResult result{};
        result.geometry.area_ratio = ar;
        result.performance.exit_mach = Me;
        result.performance.exit_velocity = ve;
        result.performance.specific_impulse = isp;
        result.performance.thrust_coefficient = 0.0; 
        result.performance.constraint_active = false;

        return result;
    }

} // namespace nozzle::design