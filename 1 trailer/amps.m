function [a1, a2, a12, a112] = amps(zf, alpha, gamma_val, h_fun)
%AMPS  The four amplitudes, a(z) = -gamma*h(z)*G^{-1}(z)*z. ONE-TRAILER
%   (r=3) version. Standalone copy of the local function inside
%   calculations.m, so other scripts (compare_gain.m, convergence_map.m)
%   can call it too -- calculations.m's own local copy still takes
%   priority when calculations.m itself runs.
    h    = h_fun(zf);
    a1   = -gamma_val*h*zf(1);
    a2   = -gamma_val*h*zf(2);
    a12  = -gamma_val*h*(zf(1)*zf(2) + zf(3));
    a112 = -gamma_val*h*(alpha*zf(2)^3 + zf(1)*zf(3) + zf(4));
end
