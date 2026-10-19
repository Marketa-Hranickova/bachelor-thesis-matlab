% COMPARE_PHYSICAL
% Same controller, two plants, ONE TRAILER (r = 3):
%   (a) the nilpotent approximation  z' = u1*g1(z) + u2*g2(z)
%       -- what simulate_case.m / calculations.m show
%   (b) the exact kinematics         x' = u1*f1(x) + u2*f2(x)
%       -- what the real robot, and Gazebo, are
% The gap is the remainder the approximation drops. Run this BEFORE
% trusting any retune and before porting anything to Gazebo: at the
% two-trailer case a setting that looked perfect on (a) diverged on (b).

P = params();
d0 = P.d0;  l1 = P.l1;  alpha = P.alpha;
gamma = P.gamma;  epsilon = P.epsilon;  omega = P.omega;
k12 = P.k12;  k112 = P.k112;
z0 = P.z0;  t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;

h_fun  = P.h_fun;                       % Lyapunov function lives in params.m

% CAREFUL: P.x0 stores the angles in DEGREES (x0_to_z0 converts them), but
% the physical simulation integrates radians.
x0 = [P.x0(1); P.x0(2); deg2rad(P.x0(3)); deg2rad(P.x0(4))];

fprintf('========================================\n');
fprintf('  x0      = [%.3f; %.3f; %.1f deg; %.1f deg]\n', ...
        P.x0(1), P.x0(2), P.x0(3), P.x0(4));
fprintf('  z0      = [%.4f; %.4f; %.4f; %.4f],  ||z0|| = %.4f\n', z0, norm(z0));
fprintf('  gamma   = %.4e,  eps = %.2f,  k = (%d,%d)\n', gamma, epsilon, k12, k112);
fprintf('  V       : %s\n', P.lyap_name);
fprintf('  time    : [%.0f, %.0f]\n', t_start, t_end);
fprintf('========================================\n');

fprintf('\nrunning approximate (z) model ...\n');
[t_a, z_a, u1_a, u2_a] = simulate_case( ...
    z0, alpha, gamma, epsilon, omega, k12, k112, ...
    t_start, t_end, dt, h_fun);

fprintf('running exact (x) model ...\n');
[t_p, x_p, z_p, u1_p, u2_p] = simulate_case_physical( ...
    x0, d0, l1, alpha, gamma, epsilon, omega, k12, k112, ...
    t_start, t_end, dt, h_fun);

% ---------------------------------------------------------------- results
na = sqrt(sum(z_a.^2, 2));   np_ = sqrt(sum(z_p.^2, 2));
hit = @(t, n, f) deal_hit(t, n, f*norm(z0));

fprintf('\n--- convergence ---\n');
for f = [0.1 0.01 0.001]
    fprintf('  ||z|| -> %6.1f%% of ||z0||   approx %10s   exact %10s\n', ...
            f*100, hit(t_a, na, f), hit(t_p, np_, f));
end
fprintf('  final ||z||                  approx %10.3e   exact %10.3e\n', ...
        na(end), np_(end));

fprintf('\n--- final physical state (exact model) ---\n');
fprintf('  pos (%+.4f, %+.4f) m   phi %+.2f deg   theta1 %+.2f deg\n', ...
        x_p(end,1), x_p(end,2), rad2deg(x_p(end,3)), rad2deg(x_p(end,4)));

fprintf('\n--- controls (same formula in both runs) ---\n');
fprintf('  approx : max|u1| = %.4f,  max|u2| = %.4f\n', max(abs(u1_a)), max(abs(u2_a)));
fprintf('  exact  : max|u1| = %.4f,  max|u2| = %.4f\n', max(abs(u1_p)), max(abs(u2_p)));
fprintf('  limits : u1_max  = %.4f,  u2_max  = %.4f\n', P.u1_max, P.u2_max);
if max(abs(u1_p)) > P.u1_max || max(abs(u2_p)) > P.u2_max
    fprintf(2, '  *** the EXACT model exceeds a limit -- do not trust this tuning ***\n');
end

% drift between the two models, on a common grid
t_c  = linspace(t_start, min(t_a(end), t_p(end)), 2000).';
drift = sqrt(sum((interp1(t_a, z_a, t_c) - interp1(t_p, z_p, t_c)).^2, 2));
fprintf('\n--- cost of the nilpotent approximation ---\n');
fprintf('  max ||z_exact - z_approx|| = %.5f  (%.1f%% of ||z0||), at t = %.1f s\n', ...
        max(drift), 100*max(drift)/norm(z0), t_c(find(drift == max(drift), 1)));

% ------------------------------------------------------------------ plots
figure('Name', 'one trailer: approximate vs exact', 'Color', 'w');

subplot(2,2,1);
semilogy(t_a, na, 'LineWidth', 1.1); hold on
semilogy(t_p, np_, 'LineWidth', 1.1);
grid on; xlabel('t [s]'); ylabel('||z||');
legend('approximate (z model)', 'exact (x model)', 'Location', 'best');
title('norm in privileged coordinates');

subplot(2,2,2);
plot(x_p(:,1), x_p(:,2), 'LineWidth', 1.1); hold on
plot(x_p(1,1), x_p(1,2), 'o'); plot(0, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);
grid on; axis equal; xlabel('x_1 [m]'); ylabel('x_2 [m]');
legend('path', 'start', 'goal', 'Location', 'best');
title('where the vehicle actually drives');

subplot(2,2,3);
plot(t_p, rad2deg(x_p(:,3:4)), 'LineWidth', 1.1);
grid on; xlabel('t [s]'); ylabel('angle [deg]');
legend('\phi', '\theta_1', 'Location', 'best');
title('physical angles, exact model');

subplot(2,2,4);
plot(t_p, u1_p, 'LineWidth', 0.8); hold on
yline( P.u1_max, '--'); yline(-P.u1_max, '--');
grid on; xlabel('t [s]'); ylabel('u_1 [m/s]');
title('driving velocity vs limit');

fprintf('\n=== Done. ===\n');


function s = deal_hit(t, n, target)
    i = find(n <= target, 1);
    if isempty(i), s = 'never'; else, s = sprintf('%.0f s', t(i)); end
end
