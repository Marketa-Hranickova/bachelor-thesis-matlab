function g = gamma_max(z, alpha, omega, k12, k112, u1max, u2max, h_fun)
%GAMMA_MAX  Largest gain for which BOTH velocity limits still hold at state z.
%   ONE-TRAILER (r=3) version. Bisection in log space: 200 halvings
%   between 1e-16 and 10.

    if nargin < 8 || isempty(h_fun)
        h_fun = @(zz) 1;            % plain quadratic Lyapunov function
    end

    lo = 1e-16;
    hi = 10;
    for i = 1:200
        mid = sqrt(lo*hi);
        [b1, b2] = control_bound(z, alpha, mid, omega, k12, k112, h_fun);
        if b1 <= u1max && b2 <= u2max
            lo = mid;
        else
            hi = mid;
        end
    end
    g = lo;
end
