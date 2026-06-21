// Copyright Roman Dial 2026
#pragma once

#include <cmath>

namespace numerics {

    /* Safeguarded Newton-Raphson root finder.
       Combines Newton steps with bisection fallback

       f, df: the function and its derivative 
       lo, hi: bracket containing the root (f must change sign across it)
       tol: stop when |f(x)| < tol
       max_iter: safety cap against non-convergence

       returns: x such that f(x) is aproxamately 0
    */
    template <typename F, typename DF>
    double safeguarded_newton(F f, DF df, double lo, double hi, double tol = 1e-10, int max_iter = 100) {
        double  f_lo = f(lo);
        double f_hi = f(hi);
        double x = 0.5 * (lo + hi);

        for (int i = 0; i < max_iter; ++i) {
            const double fx = f(x);
            const double dfx = df(x);

            if (std::abs(fx) < tol) {
                return x;
            }

            double x_new;
            const double x_newton = x - fx / dfx;

            if (dfx == 0 || x_newton < lo || x_newton > hi) {
                x_new = 0.5 * (lo + hi); // bisection fallback
            } else {
                x_new = x_newton; // newton step
            }

            const double f_new = f(x_new);

            if (f_new * f_lo > 0.0) {
                lo = x_new;
                f_lo = f_new;
            } else {
                hi = x_new;
                f_hi = f_new;
            }

            x = x_new;
        }

        return x;
    }

} // namespace numerics