function [b1, b2, terms1, terms2] = control_bound(z, alpha, kappa, theta, ...
                                                  gamma_val, omega, k12, k112, k1112, ...
                                                  h_fun)
%CONTROL_BOUND  Worst-case |u1|, |u2| at state z: the constant part plus the
%   magnitude of every oscillating term.
%
%   h_fun is the state-dependent factor of the Lyapunov function (see
%   params.m). It scales ALL five amplitudes, so leaving it out makes this
%   bound wrong by exactly that factor. Omit it only for the plain
%   quadratic V, where h = 1.


    if nargin < 10 || isempty(h_fun)
        h_fun = @(zz) 1;            % plain quadratic Lyapunov function
    end

    [a1, a2, a12, a112, a1112] = amps(z, alpha, kappa, theta, gamma_val, h_fun);

    T12   = (2*k12*omega*abs(a12))^(1/2);
    T112  = (2*k112*omega)^(2/3) * abs(a112)^(1/3);
    T1112 = (k1112*omega)^(3/4) * abs(a1112)^(1/4);

    terms1 = [abs(a1), T12,   T112, 2*T1112];
    terms2 = [abs(a2), T12, 2*T112, 6*T1112];

    b1 = sum(terms1);
    b2 = sum(terms2);
end
