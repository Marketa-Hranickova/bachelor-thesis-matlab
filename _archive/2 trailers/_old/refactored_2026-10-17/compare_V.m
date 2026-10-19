% COMPARE_V  Candidate Lyapunov functions, ranked under the robot's limits.
%
%   For a radially symmetric V = phi(r), r = ||z||, the whole function
%   enters the controller through ONE scalar factor (see params.m):
%
%       a(z) = -gamma*h(z)*G^-1(z)*z ,      h(z) = phi'(r)/(gamma*r)
%
%   so the candidates below differ ONLY in h (lyapunov_h.m). Two
%   consequences are used here:
%
%     1. Every radial V drives the state along the SAME straight ray to the
%        origin in z. Only the speed along it changes. So the admissible
%        gain can be found by walking that ray from ||z0|| down to 0.1 % of it.
%     2. The closed loop approximates the gradient flow  dr/dt = -gamma*h(r)*r ,
%        a 1-D ODE. Integrating it gives the convergence time without
%        simulating hours of oscillation. (Checked against the full
%        oscillating simulation over 80 s: they agree to 0.2 % of the
%        progress made.)
%
%   Method per candidate: the largest gamma such that
%       (1) the control bound holds at EVERY point of the ray, and
%       (2) epsilon*gamma*h stays below P.cap_product on the ray,
%   then integrate 1/(gamma*h(r)*r) for the time to reach 10 %, 1 %, 0.1 %.
%
%   The logarithmic V with c = 10 ... 100 is the one Grushkovskaya sent
%   (17 Oct 2026, trailers2.m): V = gamma/(2c)*log(1 + c*r^2).
%
%   NOTE: this uses the approximate model and the gradient-flow estimate --
%   treat the times as a ranking, then confirm the winner with
%   compare_physical.m.

P = params();
print_setup(P);

r0   = norm(P.z0);
zhat = P.z0 / r0;
rs   = r0 * logspace(-3, 0, 200);       % the ray: 0.1 % of ||z0|| up to ||z0||

% candidates: {kind, parameter}
CAND = { 'quadratic', [];
         'log',       1;
         'log',       10;
         'log',       100;
         'log',       1000;
         'sqrt',      0.005 };

fprintf('\nrays from ||z0|| = %.5f,  eps = %.1f s,  k = (%d,%d,%d)\n\n', ...
        r0, P.epsilon, P.k12, P.k112, P.k1112);
fprintf('%-24s %11s %10s %9s %9s %9s %9s\n', ...
        'Lyapunov function', 'gamma_max', 'eps*g*h', 't->10%', 't->1%', 't->0.1%', 'speedup');

base = NaN;
for c = 1:size(CAND, 1)
    [h, name] = lyapunov_h(CAND{c, 1}, CAND{c, 2});
    Q = P;
    Q.h_fun = h;

    % largest gain admissible at EVERY point of the ray (velocity limits) ...
    g   = inf;
    hmx = 0;
    for r = rs
        g   = min(g, gamma_max(r*zhat, Q));
        hmx = max(hmx, h(r*zhat));
    end
    % ... and with the contraction cap  eps*gamma*h <= cap_product
    g = min(g, P.cap_product / (P.epsilon*hmx));

    % dr/dt = -g*h(r)*r  ->  time from r0 down to frac*r0
    hr = @(r) h(r*zhat);
    T  = @(frac) integral(@(r) 1./(g*arrayfun(hr, r).*r), frac*r0, r0) / 3600;

    t10 = T(0.10);  t1 = T(0.01);  t01 = T(0.001);
    if isnan(base), base = t1; end

    fprintf('%-24s %11.3e %10.3g %7.2f h %7.2f h %7.2f h %8.2fx\n', ...
            name, g, P.epsilon*g*hmx, t10, t1, t01, base/t1);
end

fprintf(['\nNote: u1 saturates first -- the driving-velocity limit is what\n' ...
         'throttles all candidates, u2 stays well below its limit. Shifting\n' ...
         'demand from u1 to u2 (different k, see search_k.m) is a separate\n' ...
         'lever from the choice of V; in practice it is the bigger one.\n']);
