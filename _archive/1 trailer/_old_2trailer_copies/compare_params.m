% compare_params.m -- how does the control demand depend on the parameters?
%
% No integration: everything here is a substitution into the amplitude
% formula of eq. (16), so it runs instantly. Four comparisons:
%
%   1. geometry        -- the paper's two sets vs the measured hardware
%   2. period epsilon  -- demand and admissible gain against epsilon
%   3. multipliers k   -- the paper's (8,6,5) vs alternatives
%   4. term breakdown  -- which part of the control law dominates
%
% Run:  compare_params

clc;
close all;

P = params();
d0 = P.d0;  d1 = P.d1;  l1 = P.l1;  l2 = P.l2;
alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;  beta = P.beta;
gamma = P.gamma;  epsilon = P.epsilon;  omega = P.omega;
k12 = P.k12;  k112 = P.k112;  k1112 = P.k1112;
z0 = P.z0;  t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;
U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

fprintf('geometry : d0=%.3f  d1=%.3f  l1=%.3f  l2=%.3f  [m]\n', d0, d1, l1, l2);
fprintf('coeffs   : alpha=%.6f  kappa=%.6f  theta=%.6f\n', alpha, kappa, theta);


col_a = [0.00 0.20 0.70];
col_b = [0.85 0.33 0.10];
col_c = [0.00 0.55 0.30];

fprintf('\nInitial state z0 = [%g %g %g %g %g],  gamma = %g\n', z0, gamma);

% =========================================================================
%  1. GEOMETRY
% =========================================================================
geoms = { 'paper set (26)',  [0.1  0.1  1.0  1.0]
          'paper set (27)',  [1.0  1.0  1.0  1.0]
          'measured robot',  [d0   d1   l1   l2 ] };

fprintf('\n===== 1. Geometry, at eps = %.2f s, k = (%d,%d,%d) =====\n', ...
        epsilon, k12, k112, k1112);
fprintf('%-16s %9s %9s %9s %12s %12s %10s\n', ...
        'geometry', 'alpha', 'kappa', 'theta', 'max|u1|', 'max|u2|', 'gamma_max');

for i = 1:size(geoms, 1)
    g = geoms{i, 2};
    [al, ka, th] = offhooked_coeffs(g(1), g(2), g(3), g(4));
    [b1, b2] = control_bound(z0, al, ka, th, gamma, omega, k12, k112, k1112);
    gm = gamma_max(z0, al, ka, th, omega, k12, k112, k1112, U1_MAX, U2_MAX);
    fprintf('%-16s %9.5f %9.5f %9.5f %12.1f %12.1f %10.2e\n', ...
            geoms{i,1}, al, ka, th, b1, b2, gm);
end
fprintf('(max|u1| in m/s, limit %.2f;  max|u2| in rad/s, limit %.2f)\n', U1_MAX, U2_MAX);

% =========================================================================
%  2. OSCILLATION PERIOD
% =========================================================================
eps_list = [0.05 0.1 0.2 0.5 1 2 4 8 16 32];
n = numel(eps_list);
b1_eps = zeros(1, n);
b2_eps = zeros(1, n);
gm_eps = zeros(1, n);

fprintf('\n===== 2. Oscillation period (measured geometry) =====\n');
fprintf('%10s %12s %12s %12s %12s\n', 'eps [s]', 'omega', 'max|u1|', 'max|u2|', 'gamma_max');
for i = 1:n
    w = 2*pi/eps_list(i);
    [b1_eps(i), b2_eps(i)] = control_bound(z0, alpha, kappa, theta, gamma, ...
                                           w, k12, k112, k1112);
    gm_eps(i) = gamma_max(z0, alpha, kappa, theta, w, k12, k112, k1112, ...
                          U1_MAX, U2_MAX);
    fprintf('%10.2f %12.3f %12.1f %12.1f %12.3e\n', ...
            eps_list(i), w, b1_eps(i), b2_eps(i), gm_eps(i));
end

figure('Name', 'Demand and admissible gain vs period', 'NumberTitle', 'off');
subplot(1,2,1)
loglog(eps_list, b1_eps, 'o-', 'Color', col_a, 'LineWidth', 1.5); hold on
loglog(eps_list, b2_eps, 's--', 'Color', col_b, 'LineWidth', 1.5);
loglog(eps_list, U1_MAX*ones(1,n), ':', 'Color', col_a, 'LineWidth', 1.0);
loglog(eps_list, U2_MAX*ones(1,n), ':', 'Color', col_b, 'LineWidth', 1.0);
grid on; xlabel('\epsilon  [s]'); ylabel('worst-case control');
legend('|u_1| [m/s]', '|u_2| [rad/s]', 'limit u_1', 'limit u_2', 'Location', 'best');
title(sprintf('Control demand at \\gamma = %g', gamma));

subplot(1,2,2)
loglog(eps_list, gm_eps, 'o-', 'Color', col_c, 'LineWidth', 1.5);
grid on; xlabel('\epsilon  [s]'); ylabel('admissible \gamma');
title('Largest gain that fits the robot');

% =========================================================================
%  3. FREQUENCY MULTIPLIERS
% =========================================================================
% Only triples satisfying condition C1 of Theorem 1 are legitimate.
ks = [ 8  6  5      % the paper's choice
      14  5  1
      15  2  3
      12  8  1
       2  3  1 ];

fprintf('\n===== 3. Frequency multipliers (measured geometry, eps = %.2f s) =====\n', epsilon);
fprintf('%6s %6s %8s %12s %12s %12s   %s\n', ...
        'k12', 'k112', 'k1112', 'max|u1|', 'max|u2|', 'gamma_max', 'condition C1');
for i = 1:size(ks, 1)
    [b1, b2] = control_bound(z0, alpha, kappa, theta, gamma, omega, ...
                             ks(i,1), ks(i,2), ks(i,3));
    gm = gamma_max(z0, alpha, kappa, theta, omega, ks(i,1), ks(i,2), ks(i,3), ...
                   U1_MAX, U2_MAX);
    if satisfies_c1(ks(i,1), ks(i,2), ks(i,3))
        verdict = 'OK';
    else
        verdict = '*** VIOLATES C1 -- not usable ***';
    end
    fprintf('%6d %6d %8d %12.1f %12.1f %12.3e   %s\n', ...
            ks(i,1), ks(i,2), ks(i,3), b1, b2, gm, verdict);
end
fprintf(['NOTE: (2,3,1) is in the list on purpose. It gives the lowest demand\n' ...
         '      of all, and it is forbidden: 2 - 3 + 1 = 0 is a resonance of\n' ...
         '      order 3. Always run satisfies_c1 before adopting a triple.\n']);

% =========================================================================
%  4. TERM BREAKDOWN -- where does the control effort go?
% =========================================================================
[b1, b2, t1, t2] = control_bound(z0, alpha, kappa, theta, gamma, omega, ...
                                 k12, k112, k1112);

fprintf('\n===== 4. Where the control effort goes (measured geometry) =====\n');
names = {'|a|  (useful)', 'bracket 2', 'bracket 3', 'bracket 4'};
fprintf('%-16s %12s %8s %12s %8s\n', 'term', 'to |u1|', '%', 'to |u2|', '%');
for i = 1:4
    fprintf('%-16s %12.2f %7.1f%% %12.2f %7.1f%%\n', ...
            names{i}, t1(i), 100*t1(i)/b1, t2(i), 100*t2(i)/b2);
end
fprintf('%-16s %12.2f %7.1f%% %12.2f %7.1f%%\n', 'TOTAL', b1, 100, b2, 100);

figure('Name', 'Term breakdown', 'NumberTitle', 'off');
bar([t1; t2].');
set(gca, 'XTickLabel', names);
grid on; ylabel('contribution');
legend('to |u_1|', 'to |u_2|', 'Location', 'best');
title('Where the commanded velocity goes');

fprintf('\n=== Done. ===\n');
