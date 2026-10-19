function u1_val = u1_ctrl(t, zf, omega, k12, k112, alpha, gamma_val, h_fun)
%U1_CTRL  Driving velocity u1 [m/s] at time t, for the state zf frozen at the
%   start of the current sampling interval. ONE-TRAILER (r=3): two
%   oscillating terms. Amplitudes and term sizes come from osc_terms.m.
    [a1, ~, ~, ~, T12, T112] = osc_terms(zf, alpha, gamma_val, omega, k12, k112, h_fun);
    u1_val = a1 ...
        + T12  * cos(k12*omega*t) ...
        + T112 * cos(k112*omega*t);
end
