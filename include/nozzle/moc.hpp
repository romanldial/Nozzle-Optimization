// Copyright Roman Dial 2026
#pragma once

#include <vector>
#include "nozzle/geometry.hpp"  
#include "nozzle/isentropic.hpp"

namespace nozzle::moc {

    class MeshPoint {
    public:
        MeshPoint(geometry::Point2D position,
                  double theta,
                  double nu,
                  double gamma)
            : pos_{position}
            , theta_{theta}
            , nu_{nu}
            , mach_{core::mach_from_prandtl_meyer(nu, gamma)}
            , mu_{core::mach_angle(mach_)}
        {}

        const geometry::Point2D& position() const { return pos_; }
        double theta() const { return theta_; }
        double nu()    const { return nu_; }
        double mach()  const { return mach_; }
        double mu()    const { return mu_; }

    private:
        geometry::Point2D pos_;
        double theta_;
        double nu_;
        double mach_;
        double mu_;
    };

    class Solution {
        public:
        explicit Solution(double gamma)
            : gamma_{gamma}
        {}
        
        size_t add_point(MeshPoint point) {
            points_.push_back(point);
            return points_.size() - 1;
        }

        

        private:
            std::vector<MeshPoint> points_; 
            std::vector<std::vector<size_t>> characteristics_; 
            double gamma_;           
};

}