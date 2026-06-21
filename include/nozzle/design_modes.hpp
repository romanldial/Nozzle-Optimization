// Copyright Roman Dial 2026

#pragma once

#include <numbers>
#include <cmath>
#include "chamber_conditions.hpp"
#include "design_result.hpp"
#include "isentropic.hpp"

namespace nozzle::design {

    /* Design Modes for the Nozzle Optimization 
        
        Target: pe = pambient 
        
    */

    inline DesignResult design_optimal_expansion(const ChamberConditions& chamber, const double p_ambient) {
        const double p0_over_pe = chamber.p0 / p_ambient;
        const double pe_over_p0 = p_ambient / chamber.p0;
        const double Me = core::exit_mach_from_pressure_ratio(p0_over_pe, chamber.gamma);
        const double ar = core::area_ratio_from_mach(Me, chamber.gamma);
        const double ve = core::exit_velocity(chamber.T0, Me, chamber.gamma, chamber.R);
        const double isp = core::specific_impulse(ve);
        const double cf = core::thrust_coefficient(chamber.gamma, pe_over_p0, pe_over_p0, ar);
        
        DesignResult result{};
        result.geometry.area_ratio = ar;
        result.performance.exit_mach = Me;
        result.performance.exit_velocity = ve;
        result.performance.specific_impulse = isp;
        result.performance.thrust_coefficient = cf; 
        result.performance.constraint_active = false;

        return result;
    }

    inline ScaledDesignResult scale_from_exit_diameter(const DesignResult& design, const ChamberConditions& chamber, double exit_diameter) {
        const double gamma = chamber.gamma;
        const double p0 = chamber.p0;
        const double T0 = chamber.T0;
        const double R = chamber.R;

        const double exit_area = std::numbers::pi / 4.0 * exit_diameter * exit_diameter;
        const double throat_area = exit_area / design.geometry.area_ratio;
        const double throat_diameter = std::sqrt(4.0 * throat_area / std::numbers::pi);
        const double mass_flow_gamma_term = std::pow((gamma + 1.0) / 2.0, -(gamma + 1.0) / (2.0 * (gamma - 1.0)));
        const double mass_flow = (p0 * throat_area) / std::sqrt(T0) * std::sqrt(gamma / R) * mass_flow_gamma_term;
        const double thrust = design.performance.thrust_coefficient * p0 * throat_area;

        return ScaledDesignResult{
            design,
            throat_area,
            exit_area,
            throat_diameter,
            exit_diameter,
            mass_flow,
            thrust
        };

    }

    inline ScaledDesignResult design_diameter_limited(const ChamberConditions& chamber, double p_ambient, double max_exit_diameter) {
        DesignResult design = design_optimal_expansion(chamber, p_ambient);
        ScaledDesignResult scaled = scale_from_exit_diameter(design, chamber, max_exit_diameter);
        scaled.design.performance.constraint_active = false;
        return scaled;
    }

    
} // namespace nozzle::design