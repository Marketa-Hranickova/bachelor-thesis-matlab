function [g, hmax] = gamma_ray(z0, alpha, omega, epsilon, k12, k112, ...
                               u1max, u2max, h_fun, cap_product)
%GAMMA_RAY  Largest fixed gain that is admissible along the WHOLE straight
%   ray from z0 to the origin (every radial V drives the state along it).
%
%   [g, hmax] = gamma_ray(P.z0, P.alpha, P.omega, P.epsilon, P.k12, P.k112, ...
%                         P.u1_max, P.u2_max, P.h_fun, P.cap_product)
%
%   Two constraints, both must hold at every point of the ray:
%     (1) amplitude:   control_bound <= (u1max, u2max)      [gamma_max.m]
%     (2) contraction: epsilon*gamma*h(z) <= cap_product
%   g is the largest gamma meeting both; hmax is the largest h on the ray.
%   This is how the gain in params.m was chosen. Looking at z0 alone is not
%   enough for the sqrt / log V: h grows towards the origin, so (2) can
%   bind near the goal even when (1) is slack at the start.

    r0   = norm(z0);
    zhat = z0(:) / r0;
    rs   = r0 * logspace(-3, 0, 200);

    g = inf;  hmax = 0;
    for r = rs
        g    = min(g, gamma_max(r*zhat, alpha, omega, k12, k112, u1max, u2max, h_fun));
        hmax = max(hmax, h_fun(r*zhat));
    end
    g = min(g, cap_product/(epsilon*hmax));
end
