function [b1, b2, terms1, terms2] = control_bound(z, alpha, kappa, theta, ...
                                                  gamma_val, omega, k12, k112, k1112)
%CONTROL_BOUND  Worst-case magnitude of the controls (16).
%
%   Sums the absolute values of all four contributions, i.e. assumes every
%   oscillation peaks at the same instant. The true maximum over one period
%   is somewhat lower, so this is a conservative bound.
%
%   terms1 = [|a1|, T12, T112, 2*T1112]   contributions to |u1|
%   terms2 = [|a2|, T12, 2*T112, 6*T1112] contributions to |u2|

    [a1, a2, a12, a112, a1112] = controller_amps(z, alpha, kappa, theta, gamma_val);

    T12   = (2*k12*omega*abs(a12))^(1/2);
    T112  = (2*k112*omega)^(2/3) * abs(a112)^(1/3);
    T1112 = (k1112*omega)^(3/4) * abs(a1112)^(1/4);

    terms1 = [abs(a1), T12,   T112, 2*T1112];
    terms2 = [abs(a2), T12, 2*T112, 6*T1112];

    b1 = sum(terms1);
    b2 = sum(terms2);
end
