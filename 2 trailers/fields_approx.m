function [g1, g2] = fields_approx(z, alpha, kappa, theta)
%FIELDS_APPROX  Vector fields of the nilpotent approximation,
%   z' = u1*g1(z) + u2*g2(z).
%
%   The ONE place these formulas live (calculations.m, simulate_case.m,
%   simulate_case_adaptive.m and the compare_* scripts all use this file).
%   derive_symbolic.m writes the same fields symbolically.

    g1 = [1; 0; -z(2); -z(3); kappa*z(2)^3 - z(4)];
    g2 = [0; 1; 0; -alpha*z(2)^2; -theta*z(2)*z(3)];
end
