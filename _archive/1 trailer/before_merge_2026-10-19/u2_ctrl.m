function u2_val = u2_ctrl(t, zf, omega, k12, k112, alpha, gamma_val, h_fun)
%U2_CTRL  Steering velocity u2 [rad/s] at time t, for the state zf frozen at
%   the start of the current sampling interval. ONE-TRAILER (r=3): two
%   oscillating terms. Amplitudes and term sizes come from osc_terms.m.
    [~, a2, s12, s112, T12, T112] = osc_terms(zf, alpha, gamma_val, omega, k12, k112, h_fun);
    u2_val = a2 ...
        + T12 * s12 * sin(k12*omega*t) ...
        - 2*T112 * s112 * cos(2*k112*omega*t);
end
