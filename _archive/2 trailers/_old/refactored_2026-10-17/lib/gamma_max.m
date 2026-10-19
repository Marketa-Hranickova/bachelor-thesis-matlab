function g = gamma_max(z, P)
%GAMMA_MAX  Largest gain for which BOTH velocity limits still hold at state z.
%
%   g = gamma_max(z, P)
%
%   Uses P.u1_max, P.u2_max and everything control_bound uses (P.omega,
%   P.k12..., P.h_fun). To ask "what if" for another epsilon, k or V, pass
%   a modified copy:   g = gamma_max(z, set_epsilon(P, 8))
%
%   Bisection in log space: 200 halvings between 1e-16 and 10.

    lo = 1e-16;
    hi = 10;
    for i = 1:200
        mid = sqrt(lo*hi);
        [b1, b2] = control_bound(z, P, mid);
        if b1 <= P.u1_max && b2 <= P.u2_max
            lo = mid;
        else
            hi = mid;
        end
    end
    g = lo;
end
