function [a1, a2, s12, s112, T12, T112] = osc_terms( ...
        zf, alpha, gamma_val, omega, k12, k112, h_fun)
%OSC_TERMS  Everything the one-trailer controller needs at the frozen state
%   zf: the constant parts a1, a2, the signs s.. of the bracket amplitudes
%   and the sizes T.. of the two oscillating terms.
%
%   Shared by u1_ctrl.m / u2_ctrl.m (the signal) and control_bound.m (its
%   worst case), so the two can never disagree.

    [a1, a2, a12, a112] = amps(zf, alpha, gamma_val, h_fun);

    T12  = (2*k12*omega*abs(a12))^(1/2);
    T112 = (2*k112*omega)^(2/3) * abs(a112)^(1/3);

    s12 = sign(a12);  s112 = sign(a112);
end
