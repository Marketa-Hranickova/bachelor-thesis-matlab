function u1_val = u1_ctrl(t, zf, omega, k12, k112, alpha, gamma_val, h_fun)
%U1_CTRL  Standalone copy of the local function inside calculations.m.
%   ONE-TRAILER (r=3): only two oscillating terms, no k1112/a1112.
    [a1, ~, a12, a112] = amps(zf, alpha, gamma_val, h_fun);
    u1_val = a1 ...
        + (2*k12*omega*abs(a12))^(1/2) * cos(k12*omega*t) ...
        + (2*k112*omega)^(2/3) * abs(a112)^(1/3) * cos(k112*omega*t);
end
