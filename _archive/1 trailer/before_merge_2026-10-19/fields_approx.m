function [g1, g2] = fields_approx(z, alpha)
%FIELDS_APPROX  Vector fields of the nilpotent approximation of the car with
%   ONE off-hooked trailer (r = 3, 4 states):  z' = u1*g1(z) + u2*g2(z).
%
%   The ONE place these formulas live (calculations.m, the simulators and
%   the compare_* scripts all use this file). derive_symbolic.m writes the
%   same fields symbolically.

    g1 = [1; 0; -z(2); -z(3)];
    g2 = [0; 1; 0; -alpha*z(2)^2];
end
