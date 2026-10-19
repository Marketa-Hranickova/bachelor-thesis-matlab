% COMPARE_GAIN  Fixed vs. adaptive, hardware-aware gain, and how the
% admissible gain depends on epsilon. ONE-TRAILER (r=3, 4-state) version --
% see "2 trailers/compare_gain.m" for the two-trailer version.
%
%   fixed    : gamma = P.gamma for the whole run.
%   adaptive : gamma recomputed every interval = largest value that still
%              respects the robot's |u1|,|u2| limits (gamma_max.m),
%              additionally capped so eps*gamma*h(z) <= P.cap_product.
%
% Uses the Lyapunov function from params.m. Time-to-converge estimates use
% eta_ray.m (gradient flow along the ray, valid for ANY radial V); the old
% estimate (1-eps*gamma)^N assumed the quadratic V.
%
% Run:  compare_gain

clc; close all;

P = params();
alpha = P.alpha;
k12 = P.k12;  k112 = P.k112;
z0 = P.z0;  t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;
U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

h_fun  = P.h_fun;    % Lyapunov function lives in params.m

fprintf('geometry : d0=%.3f  l1=%.3f  [m]\n', P.d0, P.l1);
fprintf('z0 = [%g %g %g %g],  eps = %g,  k = (%d,%d)\n', z0, P.epsilon, k12, k112);
fprintf('V  : %s\n', P.lyap_name);

% ---------------- fixed gain ------------------------------------------------
fprintf('\n--- fixed gain, gamma = %g ---\n', P.gamma);
[t_f, z_f, u1_f, u2_f] = simulate_case(z0, alpha, P.gamma, ...
    P.epsilon, P.omega, k12, k112, t_start, t_end, dt, h_fun);
fprintf('  final norm : %.6f\n', norm(z_f(end,:)));
fprintf('  max |u1|   : %.3f  (limit %.2f)\n', max(abs(u1_f)), U1_MAX);
fprintf('  max |u2|   : %.3f  (limit %.2f)\n', max(abs(u2_f)), U2_MAX);

% ---------------- adaptive gain ---------------------------------------------
fprintf('\n--- adaptive gain, cap eps*gamma*h <= %g ---\n', P.cap_product);
[t_a, z_a, u1_a, u2_a, gamma_a] = simulate_case_adaptive(z0, alpha, ...
    P.epsilon, P.omega, k12, k112, t_start, t_end, dt, h_fun, [], [], ...
    U1_MAX, U2_MAX, P.cap_product);
fprintf('  final norm : %.6f\n', norm(z_a(end,:)));
fprintf('  max |u1|   : %.3f  (limit %.2f)\n', max(abs(u1_a)), U1_MAX);
fprintf('  max |u2|   : %.3f  (limit %.2f)\n', max(abs(u2_a)), U2_MAX);
fprintf('  gamma used : min %.4g   max %.4g   (fixed run used %.4g throughout)\n', ...
        min(gamma_a), max(gamma_a), P.gamma);

% ---------------- how long would the adaptive controller actually take? ----
% ||z|| -> 0 along the ray from z0 with the mean gain of this run.
gamma_bar = mean(gamma_a);
targets   = [0.5 0.1 0.01 0.001];
t_need    = eta_ray(z0, gamma_bar, h_fun, targets);

fprintf('\n--- estimated time to converge at the mean adaptive gain ---\n');
fprintf('  mean gamma used : %.4g   (gradient-flow estimate, eta_ray.m)\n', gamma_bar);
for i = 1:numel(targets)
    fprintf('  ||z|| down to %5.1f%% of ||z0|| : %8.0f s = %6.1f min\n', ...
            targets(i)*100, t_need(i), t_need(i)/60);
end

col_f = [0.00 0.20 0.70];
col_a = [0.85 0.33 0.10];

% ---------------- gamma_max and time-to-converge vs epsilon -----------------
eps_sweep    = [0.5 1 2 3 4 5 6 8 12 16 32];
n_sweep      = numel(eps_sweep);
gmax_sweep   = zeros(1, n_sweep);
min10_sweep  = nan(1, n_sweep);
min1_sweep   = nan(1, n_sweep);

fprintf('\n===== gamma_max and time-to-converge vs epsilon =====\n');
fprintf('(gamma = largest gain admissible along the whole ray z0 -> 0 for the velocity\n');
fprintf(' limits AND eps*gamma*h <= %.2f, see gamma_ray.m; time from eta_ray.m)\n', P.cap_product);
fprintf('%10s %14s %14s %14s %14s\n', 'eps [s]', 'gamma', 'eps*g*hmax', 'min to 10%', 'min to 1%');
for i = 1:n_sweep
    Q = set_epsilon(P, eps_sweep(i));                 % epsilon AND omega
    [gmax_sweep(i), hmax] = gamma_ray(z0, alpha, Q.omega, Q.epsilon, k12, k112, ...
                                      U1_MAX, U2_MAX, h_fun, P.cap_product);
    t = eta_ray(z0, gmax_sweep(i), h_fun, [0.1 0.01]) / 60;
    min10_sweep(i) = t(1);
    min1_sweep(i)  = t(2);
    fprintf('%10.2f %14.3e %14.3e %14.1f %14.1f\n', Q.epsilon, gmax_sweep(i), ...
            Q.epsilon*gmax_sweep(i)*hmax, min10_sweep(i), min1_sweep(i));
end
fprintf(['NOTE: the sweep uses the approximate model. params.m records eps = 2..7\n' ...
         '      as verified against the exact kinematics for this start state; for\n' ...
         '      larger eps (and for two trailers beyond eps ~ 4 s the approximation\n' ...
         '      broke down) check with compare_physical.m before believing it.\n']);


figure('Name', 'gamma_max and ETA vs epsilon (1 trailer)', 'NumberTitle', 'off');
subplot(1,2,1)
loglog(eps_sweep, gmax_sweep, 'o-', 'Color', col_a, 'LineWidth', 1.5);
grid on; xlabel('\epsilon  [s]'); ylabel('\gamma  (admissible along the ray)');
title('Largest admissible gain vs. period');

subplot(1,2,2)
loglog(eps_sweep, min10_sweep, 'o-',  'Color', col_a, 'LineWidth', 1.5); hold on
loglog(eps_sweep, min1_sweep,  's--', 'Color', col_a, 'LineWidth', 1.0);
grid on; xlabel('\epsilon  [s]'); ylabel('minutes to converge');
legend('to 10% of ||z_0||', 'to 1% of ||z_0||', 'Location', 'best');
title('Estimated time to converge at the admissible gain');

% ---------------- plot: ||z(t)|| --------------------------------------------
figure('Name', 'Fixed vs adaptive gain -- convergence (1 trailer)', 'NumberTitle', 'off');
semilogy(t_f, vecnorm(z_f, 2, 2), 'Color', col_f, 'LineWidth', 1.5); hold on
semilogy(t_a, vecnorm(z_a, 2, 2), 'Color', col_a, 'LineWidth', 1.5);
grid on; xlabel('t [s]'); ylabel('||z(t)||');
legend('fixed gain', 'adaptive gain', 'Location', 'best');
title('Convergence: fixed vs. adaptive gain');

% ---------------- plot: controls, both cases --------------------------------
figure('Name', 'Fixed vs adaptive gain -- controls (1 trailer)', 'NumberTitle', 'off');
subplot(2,1,1)
plot(t_f, u1_f, 'Color', col_f, 'LineWidth', 1.0); hold on
plot(t_a, u1_a, 'Color', col_a, 'LineWidth', 1.0);
yline(U1_MAX, 'k--'); yline(-U1_MAX, 'k--');
grid on; ylabel('u_1 [m/s]');
legend('fixed', 'adaptive', 'limit', 'Location', 'best');
title('Control inputs: fixed vs. adaptive gain');

subplot(2,1,2)
plot(t_f, u2_f, 'Color', col_f, 'LineWidth', 1.0); hold on
plot(t_a, u2_a, 'Color', col_a, 'LineWidth', 1.0);
yline(U2_MAX, 'k--'); yline(-U2_MAX, 'k--');
grid on; xlabel('t [s]'); ylabel('u_2 [rad/s]');
legend('fixed', 'adaptive', 'limit', 'Location', 'best');

% ---------------- plot: gamma(t) for the adaptive run -----------------------
figure('Name', 'Adaptive gain over time (1 trailer)', 'NumberTitle', 'off');
plot(t_a, gamma_a, 'Color', col_a, 'LineWidth', 1.5); hold on
yline(P.gamma, 'b--');
grid on; xlabel('t [s]'); ylabel('\gamma used');
legend('adaptive \gamma(t)', 'fixed \gamma', 'Location', 'best');
title('Gain chosen every interval (hardware-limited, capped \epsilon\gamma h <= cap)');

fprintf('\n=== Done. ===\n');
