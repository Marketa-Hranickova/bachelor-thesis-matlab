function [a1, a2, s12, s112, s1112, T12, T112, T1112] = osc_terms( ...
        zf, alpha, kappa, theta, gamma_val, omega, k12, k112, k1112, h_fun)
%OSC_TERMS  Everything the controller needs at the frozen state zf:
%   the constant parts a1, a2, the signs s.. of the bracket amplitudes and
%   the sizes T.. of the three oscillating terms.
%
%   Shared by u_ctrl.m (the signal) and control_bound.m (its worst case),
%   so the two can never disagree.

    [a1, a2, a12, a112, a1112] = amps(zf, alpha, kappa, theta, gamma_val, h_fun);

    T12   = (2*k12*omega*abs(a12))^(1/2);
    T112  = (2*k112*omega)^(2/3) * abs(a112)^(1/3);
    T1112 = (k1112*omega)^(3/4) * abs(a1112)^(1/4);

    s12 = sign(a12);  s112 = sign(a112);  s1112 = sign(a1112);
end
