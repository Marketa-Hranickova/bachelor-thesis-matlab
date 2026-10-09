
%   fixed    : gamma = P.gamma for the whole run.
%   adaptive : gamma recomputed every interval = largest value that still
%              respects the robot's |u1|,|u2| limits (gamma_max.m),
%              additionally capped so eps*gamma <= P.cap_product 

clc; close all;

P = params();
alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;
k12 = P.k12;  k112 = P.k112;  k1112 = P.k1112;
z0 = P.z0;  t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;
U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

h_fun  = @(z) 1;
g1_fun = @(z, alpha, kappa, theta) [1; 0; -z(2); -z(3); kappa*z(2)^3 - z(4)];
g2_fun = @(z, alpha, kappa, theta) [0; 1; 0; -alpha*z(2)^2; -theta*z(2)*z(3)];

fprintf('geometry : d0=%.3f  d1=%.3f  l1=%.3f  l2=%.3f  [m]\n', P.d0, P.d1, P.l1, P.l2);
fprintf('z0 = [%g %g %g %g %g],  eps = %g,  k = (%d,%d,%d)\n', z0, P.epsilon, k12, k112, k1112);

% ---------------- fixed gain ------------------------------------------------
fprintf('\n--- fixed gain, gamma = %g ---\n', P.gamma);
[t_f, z_f, u1_f, u2_f] = simulate_case(z0, alpha, kappa, theta, P.gamma, ...
    P.epsilon, P.omega, k12, k112, k1112, t_start, t_end, dt, h_fun, g1_fun, g2_fun);
fprintf('  final norm : %.6f\n', norm(z_f(end,:)));
fprintf('  max |u1|   : %.3f  (limit %.2f)\n', max(abs(u1_f)), U1_MAX);
fprintf('  max |u2|   : %.3f  (limit %.2f)\n', max(abs(u2_f)), U2_MAX);

% ---------------- adaptive gain ---------------------------------------------
fprintf('\n--- adaptive gain, cap eps*gamma <= %g ---\n', P.cap_product);
[t_a, z_a, u1_a, u2_a, gamma_a] = simulate_case_adaptive(z0, alpha, kappa, theta, ...
    P.epsilon, P.omega, k12, k112, k1112, t_start, t_end, dt, h_fun, g1_fun, g2_fun, ...
    U1_MAX, U2_MAX, P.cap_product);
fprintf('  final norm : %.6f\n', norm(z_a(end,:)));
fprintf('  max |u1|   : %.3f  (limit %.2f)\n', max(abs(u1_a)), U1_MAX);
fprintf('  max |u2|   : %.3f  (limit %.2f)\n', max(abs(u2_a)), U2_MAX);
fprintf('  gamma used : min %.4g   max %.4g   (fixed run used %.4g throughout)\n', ...
        min(gamma_a), max(gamma_a), P.gamma);

% ---------------- how long would the adaptive controller actually take? ----
gamma_bar   = mean(gamma_a);
eps_gamma   = P.epsilon * gamma_bar;
factor      = 1 - eps_gamma;

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
        fprintf('  ||z|| down to %5.0f%% of ||z0|| : %14.0f intervals -> %12.0f s = %8.1f days\n', ...
                tgt*100, N_int_needed, t_needed, t_needed/86400);
    end
end

col_f = [0.00 0.20 0.70];
col_a = [0.85 0.33 0.10];


eps_sweep    = [0.05 0.1 0.2 0.5 1 2 4 8 16 32 64 128];
n_sweep      = numel(eps_sweep);
gmax_sweep   = zeros(1, n_sweep);
days10_sweep = nan(1, n_sweep);
days1_sweep  = nan(1, n_sweep);

fprintf('\n===== gamma_max and time-to-converge vs epsilon =====\n');
fprintf('(measured geometry, z0 = [%g %g %g %g %g], k = (%d,%d,%d))\n', z0, k12, k112, k1112);
fprintf('%10s %14s %14s %16s %16s\n', 'eps [s]', 'gamma_max', 'eps*gamma_max', 'days to 10%', 'days to 1%');
for i = 1:n_sweep
    eps_i   = eps_sweep(i);
    omega_i = 2*pi/eps_i;
    gmax_sweep(i) = gamma_max(z0, alpha, kappa, theta, omega_i, k12, k112, k1112, U1_MAX, U2_MAX);

    eg = eps_i * gmax_sweep(i);
    f  = 1 - eg;
    if f > 0 && f < 1
        days10_sweep(i) = (log(0.1)  / log(f)) * eps_i / 86400;
        days1_sweep(i)  = (log(0.01) / log(f)) * eps_i / 86400;
    end
    fprintf('%10.2f %14.3e %14.3e %16.2f %16.2f\n', ...
            eps_i, gmax_sweep(i), eg, days10_sweep(i), days1_sweep(i));
end

figure('Name', 'gamma_max and ETA vs epsilon', 'NumberTitle', 'off');
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
figure('Name', 'Fixed vs adaptive gain -- convergence', 'NumberTitle', 'off');
semilogy(t_f, vecnorm(z_f, 2, 2), 'Color', col_f, 'LineWidth', 1.5); hold on
semilogy(t_a, vecnorm(z_a, 2, 2), 'Color', col_a, 'LineWidth', 1.5);
grid on; xlabel('t [s]'); ylabel('||z(t)||');
legend('fixed gain', 'adaptive gain', 'Location', 'best');
title('Convergence: fixed vs. adaptive gain');

% ---------------- plot: controls, both cases --------------------------------
figure('Name', 'Fixed vs adaptive gain -- controls', 'NumberTitle', 'off');
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
figure('Name', 'Adaptive gain over time', 'NumberTitle', 'off');
plot(t_a, gamma_a, 'Color', col_a, 'LineWidth', 1.5); hold on
yline(P.gamma, 'b--');
grid on; xlabel('t [s]'); ylabel('\gamma used');
legend('adaptive \gamma(t)', 'fixed \gamma', 'Location', 'best');
title('Gain chosen every interval (hardware-limited, capped eps\gamma <= cap)');

fprintf('\n=== Done. ===\n');
