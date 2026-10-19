% COMPARE_PHYSICAL  Same controller, two plants -- the honest preview of the lab.
%
%   (a) the nilpotent approximation   z' = u1*g1(z) + u2*g2(z)
%       -- what plots.m shows, the system the controller is designed on
%   (b) the exact kinematics          x' = u1*f1(x) + u2*f2(x)
%       -- what the real robot (and Gazebo) is
%
%   u1, u2 are computed by the SAME formula in both runs; only the plant
%   differs. The gap between the two is the remainder the approximation
%   drops. ALWAYS run this after retuning in params.m: a setting can look
%   perfect on (a) and still fail on (b) (see the epsilon note in params.m).
%
%   Takes a while: the exact model is integrated over the whole P.t_end.

P = params();
print_setup(P);

% =========================================================================
%  PART 0 -- the coordinate transform must be right before anything else
% =========================================================================
check_transform(P);

% =========================================================================
%  PART 1 -- both simulations, identical controller
% =========================================================================
fprintf('\nrunning approximate (z) model ...\n');
A = simulate(P, 'approx');
fprintf('running exact (x) model ...\n');
E = simulate(P, 'exact');

% =========================================================================
%  PART 2 -- comparison
% =========================================================================
z0 = P.z0;
x0 = P.x0;
xe = E.x;

fprintf('\n--- final state in z ---\n');
fprintf('  approximate : [%9.5f %9.5f %9.5f %9.5f %9.5f]  ||z|| = %.5f\n', ...
        A.z(end,:), norm(A.z(end,:)));
fprintf('  exact       : [%9.5f %9.5f %9.5f %9.5f %9.5f]  ||z|| = %.5f\n', ...
        E.z(end,:), norm(E.z(end,:)));
fprintf('  ||z0||      = %.5f\n', norm(z0));

fprintf('\n--- final state in x (what the robot physically does) ---\n');
fprintf('  start       : pos (%.4f, %.4f) m,  phi %.2f deg, th1 %.2f deg, th2 %.2f deg\n', ...
        x0(1), x0(2), rad2deg(x0(3)), rad2deg(x0(4)), rad2deg(x0(5)));
fprintf('  end         : pos (%.4f, %.4f) m,  phi %.2f deg, th1 %.2f deg, th2 %.2f deg\n', ...
        xe(end,1), xe(end,2), rad2deg(xe(end,3)), rad2deg(xe(end,4)), rad2deg(xe(end,5)));
fprintf('  ||x0|| = %.5f   ->   ||x_end|| = %.5f\n', norm(x0), norm(xe(end,:)));

% How far apart the two models drift, resampled onto a common grid.
t_c   = linspace(P.t_start, min(A.t(end), E.t(end)), 2000).';
za_c  = interp1(A.t, A.z, t_c);
ze_c  = interp1(E.t, E.z, t_c);
drift = sqrt(sum((za_c - ze_c).^2, 2));
[drift_max, i_max] = max(drift);
fprintf('\n--- how far the exact model drifts from the approximation ---\n');
fprintf('  max ||z_exact - z_approx|| : %.5f  (at t = %.2f s)\n', drift_max, t_c(i_max));
fprintf('  same, relative to ||z0||   : %.1f %%\n', 100*drift_max/norm(z0));

fprintf('\n--- controls (same formula in both runs) ---\n');
fprintf('  approximate : max |u1| = %.4f,  max |u2| = %.4f\n', max(abs(A.u1)), max(abs(A.u2)));
fprintf('  exact       : max |u1| = %.4f,  max |u2| = %.4f\n', max(abs(E.u1)), max(abs(E.u2)));
fprintf('  TurtleBot3  : u1_max  = %.4f,  u2_max  = %.4f\n', P.u1_max, P.u2_max);
if max(abs(E.u1)) > P.u1_max || max(abs(E.u2)) > P.u2_max
    fprintf('  *** the EXACT model exceeds a limit -- do not trust this tuning ***\n');
end

% =========================================================================
%  PART 3 -- plots
% =========================================================================
figure('Name', 'approximate vs exact model', 'Color', 'w');

subplot(2,2,1);
plot(A.t, sqrt(sum(A.z.^2, 2)), 'LineWidth', 1.1); hold on
plot(E.t, sqrt(sum(E.z.^2, 2)), 'LineWidth', 1.1);
grid on; xlabel('t [s]'); ylabel('||z||');
legend('approximate (z model)', 'exact (x model)', 'Location', 'best');
title('norm in privileged coordinates');

subplot(2,2,2);
plot(xe(:,1), xe(:,2), 'LineWidth', 1.1); hold on
plot(xe(1,1), xe(1,2), 'o', 'MarkerFaceColor', 'auto');
plot(0, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);
grid on; axis equal; xlabel('x_1 [m]'); ylabel('x_2 [m]');
legend('path', 'start', 'goal', 'Location', 'best');
title('where the vehicle actually drives');

subplot(2,2,3);
plot(E.t, rad2deg(xe(:,3:5)), 'LineWidth', 1.1);
grid on; xlabel('t [s]'); ylabel('angle [deg]');
legend('\phi', '\theta_1', '\theta_2', 'Location', 'best');
title('physical angles, exact model');

subplot(2,2,4);
plot(t_c, drift, 'LineWidth', 1.1);
grid on; xlabel('t [s]'); ylabel('||z_{exact} - z_{approx}||');
title('cost of the nilpotent approximation');

fprintf('\n=== Done. ===\n');
