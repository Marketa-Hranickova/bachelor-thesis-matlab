function [g1, g2] = fields_approx(z, P)
%FIELDS_APPROX  Vector fields of the NILPOTENT APPROXIMATION at state z.
%
%   [g1, g2] = fields_approx(z, P)
%
%   The approximate system is   z' = u1*g1(z) + u2*g2(z).
%   It is the truncation of the exact system in privileged coordinates --
%   the system the controller is designed on. (fields_exact.m is the real
%   plant.)  derive_symbolic.m builds the same fields symbolically and
%   checks controllability from them.

    g1 = [1; 0; -z(2); -z(3); P.kappa*z(2)^3 - z(4)];
    g2 = [0; 1; 0; -P.alpha*z(2)^2; -P.theta*z(2)*z(3)];
end
