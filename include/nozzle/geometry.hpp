// Copyright Roman Dial 2026
#pragma once

namespace nozzle::geometry {

    /* Geometric description of the Nozzle */
    struct NozzleGeometry {
        double area_ratio; // Ae/A*
    };

    /* A location in the nozzle plane. Pure geometry — no flow state. */
    struct Point2D {
        double x = 0.0; // axial coordinate  [m]
        double y = 0.0; // radial coordinate [m]
    };

} // namespace nozzle::geometry