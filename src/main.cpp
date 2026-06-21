// main.cpp — demo driver. Runs the design pipeline and prints results,
// with an inline exit-diameter round-trip sanity check.

#include <cmath>
#include <iostream>
#include <numbers>

#include "nozzle/chamber_conditions.hpp"
#include "nozzle/design_modes.hpp"
#include "nozzle/design_result.hpp"
#include "numerics/root_finding.hpp"
#include "nozzle/moc.hpp"

int main() {
    using namespace nozzle::design;
    using namespace nozzle::core;
    using namespace numerics;

    const ChamberConditions chamber{7.0e6, 3500.0, 1.22, 320.0};
    const double p_ambient = 101325.0;   // sea level, Pa
    const double exit_diameter = 1.0;     // m

    // --- Mode 2: optimal expansion ---
    const DesignResult design = nozzle::design::design_optimal_expansion(chamber, p_ambient);
    const ScaledDesignResult scaled =
        nozzle::design::scale_from_exit_diameter(design, chamber, exit_diameter);

    std::cout << "=== Mode 2: optimal expansion (De = 1.0 m) ===\n";
    std::cout << "exit Mach          " << design.performance.exit_mach << "\n";
    std::cout << "area ratio         " << design.geometry.area_ratio << "\n";
    std::cout << "thrust coefficient " << design.performance.thrust_coefficient << "\n";
    std::cout << "throat diameter    " << scaled.throat_diameter << " m\n";
    std::cout << "exit diameter      " << scaled.exit_diameter << " m\n";
    std::cout << "mass flow          " << scaled.mass_flow << " kg/s\n";
    std::cout << "thrust             " << scaled.thrust << " N\n";

    // --- Inline round-trip check ---
    // Recover De from the computed exit_area; must match the input.
    const double exit_diameter_out = std::sqrt(4.0 * scaled.exit_area / std::numbers::pi);
    const double error = std::abs(exit_diameter_out - exit_diameter);
    std::cout << "\nround-trip De: in=" << exit_diameter
              << " out=" << exit_diameter_out
              << " error=" << error
              << (error < 1.0e-9 ? "  [OK]\n" : "  [MISMATCH]\n");

    // --- Mode 1: diameter-limited (fill a 1 m bay) ---
    const ScaledDesignResult limited =
        nozzle::design::design_diameter_limited(chamber, p_ambient, exit_diameter);

    std::cout << "\n=== Mode 1: diameter-limited (D_max = 1.0 m) ===\n";
    std::cout << "throat diameter    " << limited.throat_diameter << " m\n";
    std::cout << "exit diameter      " << limited.exit_diameter << " m\n";
    std::cout << "thrust             " << limited.thrust << " N\n";
    std::cout << "constraint active  " << std::boolalpha
              << limited.design.performance.constraint_active << "\n";

    std::cout << prandtl_meyer(3.23, 1.22) * 180.0 / std::numbers::pi << "\n";

    // --- Root-finder check ---
    auto f  = [](double x){ return x * x - 2.0; };
    auto df = [](double x){ return 2.0 * x; };
    double r = numerics::safeguarded_newton(f, df, 0.0, 2.0);
    std::cout << "sqrt(2) approx " << r << "\n";

    double M = mach_from_prandtl_meyer(1.1847, 1.22);
    std::cout << "M from nu: " << M << "\n";

    nozzle::moc::MeshPoint p{ {0.0, 0.34}, 0.0, 1.185, 1.22 };
    std::cout << "M = " << p.mach() << ", mu = " << p.mu() << " rad\n";

    nozzle::moc::Solution sol{ 1.22 };
    size_t idx = sol.add_point(p);   // p is the MeshPoint you already made
    std::cout << "added point at index " << idx << "\n";

    return 0;  
}
