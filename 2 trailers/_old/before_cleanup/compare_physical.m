% COMPARE_PHYSICAL
% Same controller, two plants:
%   (a) the nilpotent approximation  z' = u1*g1(z) + u2*g2(z)
%       -- what simulate_case.m / calculations.m have been showing
%   (b) the exact kinematics         x' = u1*f1(x) + u2*f2(x)
%       -- what the real robot is
% The gap between them is the remainder the approximation drops, so this
% is the honest preview of the lab run.

P = params();
d0 = P.d0;  d1 = P.d1;  l1 = P.l1;  l2 = P.l2;
alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;
gamma = P.gamma;  epsilon = P.epsilon;  omega = P.omega;
k12 = P.k12;  k112 = P.k112;  k1112 = P.k1112;
x0 = P.x0;  z0 = P.z0;
t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;

h_fun  = @(z) 1;
g1_fun = @(z, alpha, kappa, theta) [1; 0; -z(2); -z(3); kappa*z(2)^3 - z(4)];
g2_fun = @(z, alpha, kappa, theta) [0; 1; 0; -alpha*z(2)^2; -theta*z(2)*z(3)];

% =========================================================================
%  PART 0 -- check the fast transform against the symbolic one
% =========================================================================
% x_to_z.m hard-codes the closed form of F(0); x0_to_z0.m re-derives it
% with the Symbolic Toolbox. They must agree, otherwise the hard-coded
% matrix has a typo.
fprintf('=== x_to_z vs x0_to_z0 (symbolic reference) ===\n');
rng(0);
err_fwd = 0;  err_rt = 0;
for i = 1:3
    xt = [0.5*randn; 0.5*randn; 0.8*randn; 0.8*randn; 0.8*randn];
    z_fast = x_to_z(xt, d0, d1, l1, l2);
    z_sym  = x0_to_z0(xt, d0, d1, l1, l2);
    err_fwd = max(err_fwd, max(abs(z_fast - z_sym)));
    err_rt  = max(err_rt,  max(abs(z_to_x(z_fast, d0, d1, l1, l2) - xt)));
end
fprintf('  max |x_to_z - x0_to_z0|      : %.3e\n', err_fwd);
fprintf('  max |z_to_x(x_to_z(x)) - x|  : %.3e\n', err_rt);
assert(err_fwd < 1e-9, 'x_to_z disagrees with x0_to_z0');
assert(err_rt  < 1e-9, 'z_to_x is not the inverse of x_to_z');
disp('   PASSED')

% =========================================================================
%  PART 1 -- both simulations, identical controller
% =========================================================================
fprintf('\n========================================\n');
fprintf('  x0      = [%.3f; %.3f; %.1f deg; %.1f deg; %.1f deg]\n', ...
        x0(1), x0(2), rad2deg(x0(3)), rad2deg(x0(4)), rad2deg(x0(5)));
fprintf('  z0      = [%.4f; %.4f; %.4f; %.4f; %.4f]\n', z0);
fprintf('  time    : [%.1f, %.1f],  epsilon = %.2f,  gamma = %.3f\n', ...
        t_start, t_end, epsilon, gamma);
fprintf('  k       = (%d, %d, %d)\n', k12, k112, k1112);
fprintf('========================================\n');

fprintf('\nrunning approximate (z) model ...\n');
[t_a, z_a, u1_a, u2_a] = simulate_case( ...
    z0, alpha, kappa, theta, gamma, ...
    epsilon, omega, k12, k112, k1112, ...
    t_start, t_end, dt, h_fun, g1_fun, g2_fun);

fprintf('running exact (x) model ...\n');
[t_p, x_p, z_p, u1_p, u2_p] = simulate_case_physical( ...
    x0, d0, d1, l1, l2, alpha, kappa, theta, gamma, ...
    epsilon, omega, k12, k112, k1112, ...
    t_start, t_end, dt, h_fun);

% =========================================================================
%  PART 2 -- comparison
% =========================================================================
fprintf('\n--- final state in z ---\n');
fprintf('  approximate : [%9.5f %9.5f %9.5f %9.5f %9.5f]  ||z|| = %.5f\n', ...
        z_a(end,:), norm(z_a(end,:)));
fprintf('  exact       : [%9.5f %9.5f %9.5f %9.5f %9.5f]  ||z|| = %.5f\n', ...
        z_p(end,:), norm(z_p(end,:)));
fprintf('  ||z0||      = %.5f\n', norm(z0));

fprintf('\n--- final state in x (what the robot physically does) ---\n');
fprintf('  start       : pos (%.4f, %.4f) m,  phi %.2f deg, th1 %.2f deg, th2 %.2f deg\n', ...
        x0(1), x0(2), rad2deg(x0(3)), rad2deg(x0(4)), rad2deg(x0(5)));
fprintf('  end         : pos (%.4f, %.4f) m,  phi %.2f deg, th1 %.2f deg, th2 %.2f deg\n', ...
        x_p(end,1), x_p(end,2), rad2deg(x_p(end,3)), ...
        rad2deg(x_p(end,4)), rad2deg(x_p(end,5)));
fprintf('  ||x0|| = %.5f   ->   ||x_end|| = %.5f\n', norm(x0), norm(x_p(end,:)));

% How far apart the two models drift, resampled onto a common grid.
t_c  = linspace(t_start, min(t_a(end), t_p(end)), 2000).';
za_c = interp1(t_a, z_a, t_c);
zp_c = interp1(t_p, z_p, t_c);
drift = sqrt(sum((za_c - zp_c).^2, 2));
fprintf('\n--- how far the exact model drifts from the approximation ---\n');
fprintf('  max ||z_exact - z_approx|| : %.5f  (at t = %.2f s)\n', ...
        max(drift), t_c(find(drift == max(drift), 1)));
fprintf('  same, relative to ||z0||   : %.1f %%\n', 100*max(drift)/norm(z0));

fprintf('\n--- controls (unchanged by the plant, same formula both runs) ---\n');
fprintf('  approximate : max |u1| = %.4f,  max |u2| = %.4f\n', ...
        max(abs(u1_a)), max(abs(u2_a)));
fprintf('  exact       : max |u1| = %.4f,  max |u2| = %.4f\n', ...
        max(abs(u1_p)), max(abs(u2_p)));
fprintf('  TurtleBot3  : u1_max  = %.4f,  u2_max  = %.4f\n', P.u1_max, P.u2_max);

% =========================================================================
%  PART 3 -- plots
% =========================================================================
figure('Name', 'approximate vs exact model', 'Color', 'w');

subplot(2,2,1);
plot(t_a, sqrt(sum(z_a.^2, 2)), 'LineWidth', 1.1); hold on
plot(t_p, sqrt(sum(z_p.^2, 2)), 'LineWidth', 1.1);
grid on; xlabel('t [s]'); ylabel('||z||');
legend('approximate (z model)', 'exact (x model)', 'Location', 'best');
title('norm in privileged coordinates');

subplot(2,2,2);
plot(x_p(:,1), x_p(:,2), 'LineWidth', 1.1); hold on
plot(x_p(1,1), x_p(1,2), 'o', 'MarkerFaceColor', 'auto');
plot(0, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);
grid on; axis equal; xlabel('x_1 [m]'); ylabel('x_2 [m]');
legend('path', 'start', 'goal', 'Location', 'best');
title('where the vehicle actually drives');

subplot(2,2,3);
plot(t_p, rad2deg(x_p(:,3:5)), 'LineWidth', 1.1);
grid on; xlabel('t [s]'); ylabel('angle [deg]');
legend('\phi', '\theta_1', '\theta_2', 'Location', 'best');
title('physical angles, exact model');

subplot(2,2,4);
plot(t_c, drift, 'LineWidth', 1.1);
grid on; xlabel('t [s]'); ylabel('||z_{exact} - z_{approx}||');
title('cost of the nilpotent approximation');

fprintf('\n=== Done. ===\n');
