// Copyright Roman Dial 2026
#pragma once

#include <vector>
#include <cmath>
#include <cassert>
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

    /* Strong wrappers, both hold a MeshPoint but are seperate 
       types to ensure they are passed in the correct order */
    class CPlusParent {
        public:
            explicit CPlusParent(MeshPoint p) : point_{p} {}
            const MeshPoint& point() const { return point_; }
        private:
            MeshPoint point_;
    };

    class CMinusParent {
        public:
            explicit CMinusParent(MeshPoint p) : point_{p} {}
            const MeshPoint& point() const { return point_; }
        private:
            MeshPoint point_;
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

    /* COMMENT HERE */
    inline MeshPoint interior_point(CPlusParent c_plus, CMinusParent c_minus, double gamma) {
        const MeshPoint& p1 = c_plus.point(); // C+ parent: carries theta - nu
        const MeshPoint& p2 = c_minus.point(); // C- parent: carries theta + nu

        // solve for flow at new point from two invariants along the characteristics:
        const double k_plus = p1.theta() - p1.nu();
        const double k_minus = p2.theta() + p2.nu();

        // derrive Mach and Mach Angle at new points:
        const double theta3 = 0.5 * (k_plus + k_minus);
        const double nu3 = 0.5 * (k_minus - k_plus);
        const double mach3 = core::mach_from_prandtl_meyer(nu3, gamma);
        const double mu3 = core::mach_angle(mach3);

        // averaged characteristic slopes:
        const double angle_plus = 0.5 * ((p1.theta() + p1.nu()) + (theta3 + mu3));
        const double angle_minus = 0.5 * ((p2.theta() - p2.nu()) + (theta3 - mu3));
        const double m_plus = std::tan(angle_plus);
        const double m_minus = std::tan(angle_minus);

        // locate new point as the intersection of two lines:
        const double denom = m_plus - m_minus;
        assert(std::abs(denom) > 1e-12); 

        const double x1 = p1.position().x;
        const double y1 = p1.position().y;
        const double x2 = p2.position().x;
        const double y2 = p2.position().y;

        const double x3 = ((y2 - y1) + (m_plus * x1 - m_minus * x2)) / denom;
        const double y3 = y1 + m_plus * (x3 - x1);

        // assemble solved node:
        return MeshPoint{ geometry::Point2D{x3, y3}, theta3, nu3, gamma };
    }

    /* COMMENT HERE */
    inline MeshPoint axis_point(CMinusParent c_minus, double gamma) {
        // flow solve
        const MeshPoint& p1 = c_minus.point(); 
        const double k_minus = p1.theta() + p1.nu(); // invariant carried down the C- line
        const double theta3 = 0.0; // boundary condition where axis is purely axial
        const double nu3 = k_minus;

        // derrive Mach and Mach Angle at new points:
        const double mach3 = core::mach_from_prandtl_meyer(nu3, gamma);
        const double mu3 = core::mach_angle(mach3);

        // averaged C- slope:
        const double angle_minus = 0.5 * ((p1.theta() - p1.mu()) + (theta3 - mu3));
        const double m_minus = std::tan(angle_minus);

        // intersect C- line with y = 0 axis:
        const double x1 = p1.position().x;
        const double y1 = p1.position().y;

        assert(std::abs(m_minus) > 1e-12);
        const double x3 = x1 - y1 / m_minus;
        const double y3 = 0.0;

        // assemble solved node:
        return MeshPoint{ geometry::Point2D{x3, y3}, theta3, nu3, gamma };
    }

    /* COMMENT HERE */
    inline MeshPoint wall_point(CPlusParent c_plus, const MeshPoint& prev_wall, double gamma) {
        const MeshPoint& p1 = c_plus.point();

        // wall segment slope:
        const double theta_wall = prev_wall.theta(); 
        const double m_wall = std::tan(theta_wall);

        // C+ characteristic slope from the parent:
        const double angle_plus = p1.theta() + p1.mu();
        const double m_plus = std::tan(angle_plus);

        // intersect C+ line with wall segment:
        const double denom = m_plus - m_wall;
        assert(std::abs(denom) > 1e-12);   // for C+ nearly parallel to wall 

        const double x1 = p1.position().x;
        const double y1 = p1.position().y;
        const double xw = prev_wall.position().x;
        const double yw = prev_wall.position().y;

        const double x3 = ((yw - y1) + (m_plus * x1 - m_wall * xw)) / denom;
        const double y3 = y1 + m_plus * (x3 - x1);

        // flow at wall point:
        const double theta3 = theta_wall; 
        const double k_plus = p1.theta() - p1.nu(); 
        const double nu3    = theta3 - k_plus; 

        // assemble wall point:
        return MeshPoint{ geometry::Point2D{x3, y3}, theta3, nu3, gamma };
    }
}