// Copyright Roman Dial 2026

#include <iostream>

#include "nozzle/chamber_conditions.hpp"
#include "nozzle/design_modes.hpp"

int main() {
    // A sample design point: a kerolox-ish chamber, sea-level target.
    const nozzle::design::ChamberConditions chamber{
        .p0    = 7.0e6,   // 7 MPa chamber pressure        [Pa]
        .T0    = 3500.0,  // 3500 K flame temperature      [K]
        .gamma = 1.22,    // ratio of specific heats        [-]
        .R     = 320.0    // specific gas constant     [J/(kg·K)]
    };

    const double p_ambient = 101325.0;  // sea level [Pa]

    const nozzle::design::DesignResult result =
        nozzle::design::design_optimal_expansion(chamber, p_ambient);

    std::cout << "--- Optimal-expansion design (sea level) ---\n";
    std::cout << "Exit Mach number     : " << result.performance.exit_mach        << "\n";
    std::cout << "Expansion ratio Ae/A*: " << result.geometry.area_ratio          << "\n";
    std::cout << "Exit velocity  [m/s] : " << result.performance.exit_velocity    << "\n";
    std::cout << "Specific impulse [s] : " << result.performance.specific_impulse << "\n";

    return 0;
}