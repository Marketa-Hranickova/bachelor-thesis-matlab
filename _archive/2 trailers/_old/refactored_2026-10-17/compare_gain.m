% COMPARE_GAIN  Fixed gain versus adaptive gain, and how epsilon limits both.
%
%   fixed    : gamma = P.gamma for the whole run.
%   adaptive : gamma recomputed every interval = the largest value that
%              still respects the robot's |u1|, |u2| limits (gamma_max.m),
%              additionally capped so that epsilon*gamma <= P.cap_product.
%
%   Parts:  1 simulate both   2 how long would adaptive take (estimate)
%           3 gamma_max versus epsilon   4 plots
%
%   NOTE: this is the APPROXIMATE model. Anything that looks good here
%   still has to survive compare_physical.m.

clc; close all;

P = params();
print_setup(P);

col_f = [0.00 0.20 0.70];       % fixed    -- blue
col_a = [0.85 0.33 0.10];       % adaptive -- orange

% =========================================================================
%  PART 1 -- the two simulations
% =========================================================================
fprintf('\n--- fixed gain, gamma = %g ---\n', P.gamma);
F = simulate(P, 'approx', 'fixed');
fprintf('  final norm : %.6f\n', norm(F.z(end,:)));
fprintf('  max |u1|   : %.3f  (limit %.2f)\n', max(abs(F.u1)), P.u1_max);
fprintf('  max |u2|   : %.3f  (limit %.2f)\n', max(abs(F.u2)), P.u2_max);

fprintf('\n--- adaptive gain, cap eps*gamma <= %g ---\n', P.cap_product);
A = simulate(P, 'approx', 'adaptive');
fprintf('  final norm : %.6f\n', norm(A.z(end,:)));
fprintf('  max |u1|   : %.3f  (limit %.2f)\n', max(abs(A.u1)), P.u1_max);
fprintf('  max |u2|   : %.3f  (limit %.2f)\n', max(abs(A.u2)), P.u2_max);
fprintf('  gamma used : min %.4g   max %.4g   (fixed run used %.4g throughout)\n', ...
        min(A.gamma), max(A.gamma), P.gamma);

% =========================================================================
%  PART 2 -- how long would the adaptive controller actually take?
% =========================================================================
% One interval shrinks ||z|| by the factor (1 - eps*gamma) (the proof's
% linear bound), so reaching a fraction f takes log(f)/log(1-eps*gamma)
% intervals.
gamma_bar = mean(A.gamma);
eps_gamma = P.epsilon * gamma_bar;
factor    = 1 - eps_gamma;

fprintf('\n--- estimated time to converge at this adaptive gain ---\n');
fprintf('  mean gamma used : %.4g   eps*gamma : %.4g   per-interval factor (1-eps*gamma) : %.6f\n', ...
        gamma_bar, eps_gamma, factor);

if ~strcmp(P.lyap_name, 'quadratic')
    fprintf(['  CAUTION: this linear estimate assumes h = 1. With V = %s the pull grows\n' ...
             '  near the goal, so the true time is much SHORTER -- compare_V.m accounts for h(z).\n'], ...
            P.lyap_name);
end

if factor <= 0
    fprintf('  eps*gamma >= 1: the linear estimate itself does not contract here -- no ETA.\n');
elseif factor >= 1
    fprintf('  factor ~ 1: convergence is immeasurably slow at this resolution.\n');
else
    for target = [0.5 0.1 0.01]
        n_needed = log(target) / log(factor);
        t_needed = n_needed * P.epsilon;
        fprintf('  ||z|| down to %5.0f%% of ||z0|| : %14.0f intervals -> %12.0f s = %8.1f days\n', ...
                target*100, n_needed, t_needed, t_needed/86400);
    end
end

% =========================================================================
%  PART 3 -- largest hardware-safe gain versus epsilon
% =========================================================================
eps_sweep = [0.05 0.1 0.2 0.5 1 2 4 8 16 32 64 128];
n_sweep   = numel(eps_sweep);
gmax_sw   = zeros(1, n_sweep);
days10_sw = nan(1, n_sweep);
days1_sw  = nan(1, n_sweep);

fprintf('\n===== gamma_max and time-to-converge vs epsilon =====\n');
fprintf('(z0 = [%g %g %g %g %g], k = (%d,%d,%d))\n', P.z0, P.k12, P.k112, P.k1112);
fprintf('%10s %14s %14s %16s %16s\n', 'eps [s]', 'gamma_max', 'eps*gamma_max', 'days to 10%', 'days to 1%');
for i = 1:n_sweep
    eps_i = eps_sweep(i);
    Q     = set_epsilon(P, eps_i);                 % keeps omega consistent
    gmax_sw(i) = gamma_max(P.z0, Q);

    eg = eps_i * gmax_sw(i);
    f  = 1 - eg;
    if f > 0 && f < 1
        days10_sw(i) = (log(0.1)  / log(f)) * eps_i / 86400;
        days1_sw(i)  = (log(0.01) / log(f)) * eps_i / 86400;
    end
    fprintf('%10.2f %14.3e %14.3e %16.2f %16.2f\n', ...
            eps_i, gmax_sw(i), eg, days10_sw(i), days1_sw(i));
end
fprintf('(large eps looks fast here -- but P.epsilon is capped by the VALIDITY of\n');
fprintf(' the approximation, see the notes in params.m and compare_physical.m;\n');
fprintf(' and the day counts assume h = 1, see the CAUTION above)\n');

% =========================================================================
%  PART 4 -- plots
% =========================================================================
figure('Name', 'gamma_max and ETA vs epsilon', 'NumberTitle', 'off');
subplot(1,2,1)
loglog(eps_sweep, gmax_sw, 'o-', 'Color', col_a, 'LineWidth', 1.5);
grid on; xlabel('\epsilon  [s]'); ylabel('\gamma_{max}  (hardware-safe)');
title('Largest hardware-safe gain vs. period');

subplot(1,2,2)
loglog(eps_sweep, days10_sw, 'o-',  'Color', col_a, 'LineWidth', 1.5); hold on
loglog(eps_sweep, days1_sw,  's--', 'Color', col_a, 'LineWidth', 1.0);
grid on; xlabel('\epsilon  [s]'); ylabel('days to converge');
legend('to 10% of ||z_0||', 'to 1% of ||z_0||', 'Location', 'best');
title('Estimated time to converge at gamma_{max}(\epsilon)');

figure('Name', 'Fixed vs adaptive gain -- convergence', 'NumberTitle', 'off');
semilogy(F.t, sqrt(sum(F.z.^2, 2)), 'Color', col_f, 'LineWidth', 1.5); hold on
semilogy(A.t, sqrt(sum(A.z.^2, 2)), 'Color', col_a, 'LineWidth', 1.5);
grid on; xlabel('t [s]'); ylabel('||z(t)||');
legend('fixed gain', 'adaptive gain', 'Location', 'best');
title('Convergence: fixed vs. adaptive gain');

figure('Name', 'Fixed vs adaptive gain -- controls', 'NumberTitle', 'off');
subplot(2,1,1)
plot(F.t, F.u1, 'Color', col_f, 'LineWidth', 1.0); hold on
plot(A.t, A.u1, 'Color', col_a, 'LineWidth', 1.0);
yline(P.u1_max, 'k--'); yline(-P.u1_max, 'k--');
grid on; ylabel('u_1 [m/s]');
legend('fixed', 'adaptive', 'limit', 'Location', 'best');
title('Control inputs: fixed vs. adaptive gain');

subplot(2,1,2)
plot(F.t, F.u2, 'Color', col_f, 'LineWidth', 1.0); hold on
plot(A.t, A.u2, 'Color', col_a, 'LineWidth', 1.0);
yline(P.u2_max, 'k--'); yline(-P.u2_max, 'k--');
grid on; xlabel('t [s]'); ylabel('u_2 [rad/s]');
legend('fixed', 'adaptive', 'limit', 'Location', 'best');

figure('Name', 'Adaptive gain over time', 'NumberTitle', 'off');
plot(A.t, A.gamma, 'Color', col_a, 'LineWidth', 1.5); hold on
yline(P.gamma, 'b--');
grid on; xlabel('t [s]'); ylabel('\gamma used');
legend('adaptive \gamma(t)', 'fixed \gamma', 'Location', 'best');
title('Gain chosen every interval (hardware-limited, capped eps\gamma <= cap)');

fprintf('\n=== Done. ===\n');
