function t = eta_ray(z0, gamma_val, h_fun, fracs)
%ETA_RAY  Time [s] until ||z|| has dropped to frac*||z0||, for any radial V.
%
%   t = eta_ray(P.z0, P.gamma, P.h_fun, [0.1 0.01 0.001])
%
%   For a radial V the closed loop follows (on average over the
%   oscillations) the gradient flow   dr/dt = -gamma*h(r)*r,   r = ||z||,
%   along the straight ray from z0 to the origin. So
%       t(frac) = integral from frac*r0 to r0 of  dr / (gamma*h(r)*r).
%   For the quadratic V (h = 1) this is ln(1/frac)/gamma; for the log and
%   sqrt V it is NOT, which is why a linear estimate 1-eps*gamma per
%   interval is wrong for them.
%
%   (Checked against the full oscillating simulation over 80 s: the two
%   agree to 0.2% of the progress made -- see compare_V.m.)

    r0   = norm(z0);
    zhat = z0(:) / r0;
    hr   = @(r) arrayfun(@(rr) h_fun(rr*zhat), r);

    t = zeros(size(fracs));
    for i = 1:numel(fracs)
        t(i) = integral(@(r) 1 ./ (gamma_val*hr(r).*r), fracs(i)*r0, r0);
    end
end
