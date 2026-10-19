function [b1, b2, terms1, terms2] = control_bound(z, alpha, gamma_val, ...
                                                  omega, k12, k112, h_fun)
%CONTROL_BOUND  Worst-case magnitude of the controls, r=3 (1 trailer).
%   Sums the absolute values of all contributions (peaks-coincide bound,
%   see the 2-trailer control_bound.m). One fewer term than the 2-trailer
%   version -- no k1112/a1112.
%
%   terms1 = [|a1|, T12, T112]     contributions to |u1|
%   terms2 = [|a2|, T12, 2*T112]   contributions to |u2|

    % h_fun is the state-dependent factor of the Lyapunov function (see
    % params.m). It scales ALL amplitudes, so leaving it out makes this
    % bound wrong by exactly that factor.
    if nargin < 7 || isempty(h_fun)
        h_fun = @(zz) 1;            % plain quadratic Lyapunov function
    end
    [a1, a2, ~, ~, T12, T112] = osc_terms(z, alpha, gamma_val, omega, k12, k112, h_fun);

    terms1 = [abs(a1), T12,   T112];
    terms2 = [abs(a2), T12, 2*T112];

    b1 = sum(terms1);
    b2 = sum(terms2);
end
