% compare_gain.m -- fixed gain (calculations.m) vs. adaptive, hardware-aware
% gain. ONE-TRAILER (r=3, 4-state) version -- see "2 trailers/compare_gain.m"
% for the two-trailer version.
%
%   fixed    : gamma = P.gamma for the whole run.
%   adaptive : gamma recomputed every interval = largest value that still
%              respects the robot's |u1|,|u2| limits (gamma_max.m),
%              additionally capped so eps*gamma <= P.cap_product.
%
% Run:  compare_gain

clc; close all;

P = params();
alpha = P.alpha;
k12 = P.k12;  k112 = P.k112;
z0 = P.z0;  t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;
U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

h_fun  = @(z) 1;
g1_fun = @(z, alpha) [1; 0; -z(2); -z(3)];
g2_fun = @(z, alpha) [0; 1; 0; -alpha*z(2)^2];

fprintf('geometry : d0=%.3f  l1=%.3f  [m]\n', P.d0, P.l1);
fprintf('z0 = [%g %g %g %g],  eps = %g,  k = (%d,%d)\n', z0, P.epsilon, k12, k112);

% ---------------- fixed gain ------------------------------------------------
fprintf('\n--- fixed gain, gamma = %g ---\n', P.gamma);
[t_f, z_f, u1_f, u2_f] = simulate_case(z0, alpha, P.gamma, ...
    P.epsilon, P.omega, k12, k112, t_start, t_end, dt, h_fun, g1_fun, g2_fun);
fprintf('  final norm : %.6f\n', norm(z_f(end,:)));
fprintf('  max |u1|   : %.3f  (limit %.2f)\n', max(abs(u1_f)), U1_MAX);
fprintf('  max |u2|   : %.3f  (limit %.2f)\n', max(abs(u2_f)), U2_MAX);

% ---------------- adaptive gain ---------------------------------------------
fprintf('\n--- adaptive gain, cap eps*gamma <= %g ---\n', P.cap_product);
[t_a, z_a, u1_a, u2_a, gamma_a] = simulate_case_adaptive(z0, alpha, ...
    P.epsilon, P.omega, k12, k112, t_start, t_end, dt, h_fun, g1_fun, g2_fun, ...
    U1_MAX, U2_MAX, P.cap_product);
fprintf('  final norm : %.6f\n', norm(z_a(end,:)));
fprintf('  max |u1|   : %.3f  (limit %.2f)\n', max(abs(u1_a)), U1_MAX);
fprintf('  max |u2|   : %.3f  (limit %.2f)\n', max(abs(u2_a)), U2_MAX);
fprintf('  gamma used : min %.4g   max %.4g   (fixed run used %.4g throughout)\n', ...
        min(gamma_a), max(gamma_a), P.gamma);

% ---------------- how long would the adaptive controller actually take? ----
gamma_bar = mean(gamma_a);
eps_gamma = P.epsilon * gamma_bar;
factor    = 1 - eps_gamma;

fprintf('\n--- estimated time to converge at this adaptive gain ---\n');
fprintf('  mean gamma used : %.4g   eps*gamma : %.4g   per-interval factor (1-eps*gamma) : %.6f\n', ...
        gamma_bar, eps_gamma, factor);

if factor <= 0
    fprintf('  eps*gamma >= 1: the linear estimate itself does not contract here -- no ETA.\n');
elseif factor >= 1
    fprintf('  factor ~ 1: convergence is immeasurably slow at this resolution.\n');
else
    targets = [0.5 0.1 0.01];
    for tgt = targets
        N_int_needed = log(tgt) / log(factor);
        t_needed     = N_int_needed * P.epsilon;
        fprintf('  ||z|| down to %5.0f%% of ||z0|| : %14.0f intervals -> %12.0f s = %8.1f minutes\n', ...
                tgt*100, N_int_needed, t_needed, t_needed/60);
    end
end

col_f = [0.00 0.20 0.70];
col_a = [0.85 0.33 0.10];


figure('Name', 'gamma_max and ETA vs epsilon (1 trailer)', 'NumberTitle', 'off');
subplot(1,2,1)
loglog(eps_sweep, gmax_sweep, 'o-', 'Color', col_a, 'LineWidth', 1.5);
grid on; xlabel('\epsilon  [s]'); ylabel('\gamma_{max}  (hardware-safe)');
title('Largest hardware-safe gain vs. period');

subplot(1,2,2)
loglog(eps_sweep, days10_sweep, 'o-',  'Color', col_a, 'LineWidth', 1.5); hold on
loglog(eps_sweep, days1_sweep,  's--', 'Color', col_a, 'LineWidth', 1.0);
grid on; xlabel('\epsilon  [s]'); ylabel('days to converge');
legend('to 10% of ||z_0||', 'to 1% of ||z_0||', 'Location', 'best');
title('Estimated time to converge at gamma_{max}(\epsilon)');

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
title('Gain chosen every interval (hardware-limited, capped eps\gamma <= cap)');

fprintf('\n=== Done. ===\n');
