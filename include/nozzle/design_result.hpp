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
    
} // namespace nozzle::design