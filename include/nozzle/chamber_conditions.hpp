#pragma once 

// Copyright Roman Dial 2026

namespace nozzle::design {

    /* Intensive Inputs that describe the chamber conditions
        
        p0: chamber pressure [Pa]
        T0: chamber temperature [K]
        gamma: specific heat ratio [-]
        R: specific gas constant [J/(kg*K)]
    
    */

struct ChamberConditions {
    double p0;
    double T0;
    double gamma;
    double R;
};

} // namespace nozzle::design