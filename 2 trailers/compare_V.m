% COMPARE_V  Candidate Lyapunov functions, ranked under the robot's limits.
%
% For a radially symmetric V = phi(r), r = ||z||, the whole function enters
% the controller through one scalar factor (see params.m):
%
%       a(z) = -gamma*h(z)*G^-1(z)*z,      h(z) = phi'(r)/(gamma*r)
%
% so the candidates below differ ONLY in h. Two consequences used here:
%
%   1. Every radial V drives the state along the SAME straight ray to the
%      origin in z; only the speed along it changes. So the admissible gain
%      can be found by walking that ray.
%   2. The closed loop approximates the gradient flow  dr/dt = -gamma*h(r)*r,
%      a 1-D ODE, so the convergence time follows without simulating hours
%      of oscillation (eta_ray.m).
%
% Per candidate, gamma is the largest value that keeps
%   (1) the control bound within the limits at EVERY point of the ray, and
%   (2) eps*gamma*h <= P.cap_product on the ray (contraction, see params.m).
% The row marked "in use" is whatever params.m currently sets.

P = params();
alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;
omega = P.omega;  k12 = P.k12;  k112 = P.k112;  k1112 = P.k1112;
U1 = P.u1_max;  U2 = P.u2_max;
z0 = P.z0;  r0 = norm(z0);  zhat = z0/r0;

% {label, h}
CAND = { ...
    ['in use: ' P.lyap_name],               P.h_fun; ...
    'quadratic  V = gamma*r^2/2',           @(z) 1; ...
    'log c=1    V = gamma/(2c)log(1+cr^2)', @(z) 1/(1 +   1*norm(z)^2); ...
    'log c=10   (Grushkovskaya''s start)',  @(z) 1/(1 +  10*norm(z)^2); ...
    'log c=100',                            @(z) 1/(1 + 100*norm(z)^2); ...
    'sqrt d=0.01 V = gamma(sqrt(r^2+d^2)-d)', @(z) 1/sqrt(norm(z)^2 + 0.01^2); ...
    'sqrt d=0.005',                         @(z) 1/sqrt(norm(z)^2 + 0.005^2) };

rs = r0 * logspace(-3, 0, 200);    % the ray, from 0.1% of ||z0|| up to it

fprintf('||z0|| = %.5f,  limits |u1|<=%.2f, |u2|<=%.2f,  k=(%d,%d,%d), eps=%.1f, cap=%.2f\n\n', ...
        r0, U1, U2, k12, k112, k1112, P.epsilon, P.cap_product);
nC  = size(CAND, 1);
res = zeros(nC, 7);     % gamma, t10, t1, t01, eps*g*hmax, max u1 use, max u2 use
for c = 1:nC
    h = CAND{c, 2};

    % (1) largest gain admissible at EVERY point of the ray
    g = inf;
    hmax = 0;
    for r = rs
        g    = min(g, gamma_max(r*zhat, alpha, kappa, theta, omega, ...
                                k12, k112, k1112, U1, U2, h));
        hmax = max(hmax, h(r*zhat));
    end
    % (2) contraction cap
    g = min(g, P.cap_product/(P.epsilon*hmax));

    % how much of each limit the chosen gain uses (worst point of the ray)
    use1 = 0;  use2 = 0;
    for r = rs
        [b1, b2] = control_bound(r*zhat, alpha, kappa, theta, g, omega, ...
                                 k12, k112, k1112, h);
        use1 = max(use1, b1/U1);  use2 = max(use2, b2/U2);
    end

    t = eta_ray(z0, g, h, [0.10 0.01 0.001]) / 3600;      % hours
    res(c, :) = [g, t, P.epsilon*g*hmax, use1, use2];
end
base = res(2, 3);                    % row 2 = quadratic V, time to 1%

fprintf('%-40s %10s %8s %8s %8s %8s %9s %6s %6s\n', 'Lyapunov function', 'gamma', ...
        't->10%', 't->1%', 't->0.1%', 'speedup', 'eps*g*h', 'u1 %', 'u2 %');
for c = 1:nC
    fprintf('%-40s %10.3e %6.2f h %6.2f h %6.2f h %7.2fx %9.2e %5.0f%% %5.0f%%\n', ...
            CAND{c,1}, res(c,1), res(c,2), res(c,3), res(c,4), base/res(c,3), ...
            res(c,5), 100*res(c,6), 100*res(c,7));
end
fprintf(['\nspeedup = time to 1%% of the quadratic V divided by this row''s.\n' ...
         'eps*g*h = largest eps*gamma*h on the ray (contraction, must stay <= cap).\n' ...
         'u1 %%, u2 %% = largest share of each velocity limit used on the ray;\n' ...
         '100%% marks the limit that sets gamma.\n']);
