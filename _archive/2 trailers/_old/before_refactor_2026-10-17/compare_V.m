% COMPARE_V
% Candidate Lyapunov functions, ranked under the robot's velocity limits.
%
% For a radially symmetric V = phi(r), r = ||z||, the whole function enters
% the controller through one scalar factor (see params.m):
%
%       a(z) = -gamma*h(z)*G^-1(z)*z,      h(z) = phi'(r)/(gamma*r)
%
% so the three candidates below differ ONLY in h. Two consequences used
% here:
%
%   1. Every radial V drives the state along the SAME straight ray to the
%      origin in z. Only the speed along it changes. So the admissible gain
%      can be found by walking that ray.
%   2. The closed loop approximates the gradient flow  dr/dt = -gamma*h(r)*r,
%      which is a 1-D ODE. Integrating it gives the convergence time without
%      simulating hours of oscillation. (Checked against the full
%      oscillating simulation over 80 s: the two agree to 0.2% of the
%      progress made.)
%
% Method per candidate: find the largest gamma whose control bound holds at
% EVERY point of the ray, then integrate 1/(gamma*h(r)*r).

P = params();
alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;
omega = P.omega;  k12 = P.k12;  k112 = P.k112;  k1112 = P.k1112;
U1 = P.u1_max;  U2 = P.u2_max;
z0 = P.z0;  r0 = norm(z0);  zhat = z0/r0;

delta = 0.01;                      % smoothing of the ||z|| candidate

CAND = { ...
    'V = gamma*r^2/2  (current)',   @(z) 1; ...
    'V = (gamma/2)*log(1+r^2)',     @(z) 1/(1 + norm(z)^2); ...
    'V = gamma*(sqrt(r^2+d^2)-d)',  @(z) 1/sqrt(norm(z)^2 + delta^2) };

rs = r0 * logspace(-3, 0, 60);     % the ray, from 0.1% of ||z0|| up to it

fprintf('||z0|| = %.5f,  limits |u1|<=%.2f, |u2|<=%.2f,  k=(%d,%d,%d), eps=%.1f\n\n', ...
        r0, U1, U2, k12, k112, k1112, P.epsilon);
fprintf('%-32s %11s %9s %9s %9s %9s\n', ...
        'Lyapunov function', 'gamma', 't->10%', 't->1%', 't->0.1%', 'speedup');

base = NaN;
for c = 1:size(CAND, 1)
    h = CAND{c, 2};

    % largest gain admissible at EVERY point of the ray
    g = inf;
    for r = rs
        g = min(g, gamma_max(r*zhat, alpha, kappa, theta, omega, ...
                             k12, k112, k1112, U1, U2, h));
    end

    % dr/dt = -g*h(r)*r  ->  time from r0 down to frac*r0
    hr = @(r) h(r*zhat);
    T  = @(frac) integral(@(r) 1./(g*arrayfun(hr, r).*r), frac*r0, r0) / 3600;

    t10 = T(0.10);  t1 = T(0.01);  t01 = T(0.001);
    if isnan(base), base = t1; end

    fprintf('%-32s %11.3e %7.2f h %7.2f h %7.2f h %8.2fx\n', ...
            CAND{c,1}, g, t10, t1, t01, base/t1);
end

fprintf(['\nNote: u1 saturates first -- at these gains |u2| peaks near 0.53 of\n' ...
         'the 2.84 rad/s available, so the driving-velocity limit is what\n' ...
         'throttles all three. Shifting demand from u1 to u2 (different k)\n' ...
         'is a separate lever from the choice of V.\n']);
