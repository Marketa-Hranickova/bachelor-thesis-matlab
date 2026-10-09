function g = gamma_max(z, alpha, kappa, theta, omega, k12, k112, k1112, u1max, u2max)
%GAMMA_MAX  Largest gain for which BOTH velocity limits still hold at state z.
%   Bisection in log space: 200 halvings between 1e-16 and 10.

    lo = 1e-16;
    hi = 10;
    for i = 1:200
        mid = sqrt(lo*hi);
        [b1, b2] = control_bound(z, alpha, kappa, theta, mid, omega, k12, k112, k1112);
        if b1 <= u1max && b2 <= u2max
            lo = mid;
        else
            hi = mid;
        end
    end
    g = lo;
end
