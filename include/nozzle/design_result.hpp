// Copyright Roman Dial 2026

#pragma once

namespace nozzle::geometry {

    /* Geometric description of the Nozzle */
    struct NozzleGeometry {
        double area_ratio; // Ae/A*
    };

} // namespace nozzle::geometry

namespace nozzle::design {

    /* Design performance of the Nozzle */
    struct NozzlePerformance {
        double exit_mach; // exit Mach number [-]
        double exit_velocity; // exit velocity [m/s]
        double specific_impulse; // specific impulse [s]
        double thrust_coefficient; // thrust coefficient [-]
        bool constraint_active; // did a constraint 
    };

    /* Complete result of a design run*/
    struct DesignResult {
        geometry::NozzleGeometry geometry;
        design::NozzlePerformance performance;
    };

    /* Complete design result with scaling */
    struct ScaledDesignResult {
        DesignResult design; 

        double throat_area; // A*  [m^2]
        double exit_area; // Ae  [m^2]
        double throat_diameter; // D*  [m]
        double exit_diameter; // De  [m]
        double mass_flow; // mdot [kg/s], ideal/choked
        double thrust; // F   [N]
    };
    
} // namespace nozzle::design