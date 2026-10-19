function [a1, a2, a12, a112] = controller_amps(z, alpha, gamma_val)
%CONTROLLER_AMPS  The four amplitudes a(z) = -gamma * G^{-1}(z) * z, h = 1.
%   ONE-TRAILER (r=3) version -- see calculations.m for the derivation and
%   the caveat about this control-law shape being an unconfirmed guess.

    a1   = -gamma_val*z(1);
    a2   = -gamma_val*z(2);
    a12  = -gamma_val*(z(1)*z(2) + z(3));
    a112 = -gamma_val*(alpha*z(2)^3 + z(1)*z(3) + z(4));
end
