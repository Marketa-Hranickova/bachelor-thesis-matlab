% plots.m -- figures for the ONE-TRAILER (r=3, 4-state) simulation.
% Run calculations.m first (it leaves t_sim, z_sim, u1_sim, u2_sim, P in
% the workspace).
%
% Run:  plots

P = params();
U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

col     = [0.00 0.20 0.70];     % approximate (nilpotent) model
col_ex  = [0.85 0.33 0.10];     % exact kinematics

% If calculations.m has not been run in this session, run it -- this used
% to be the usual reason plots.m "stopped working".
if ~exist('z_sim', 'var') || ~exist('u1_sim', 'var')
    fprintf('(t_sim/z_sim not in the workspace -- running calculations.m first)\n');
    calculations;
end

% The SAME controller on the exact kinematics. u1, u2 are not "converted"
% here: the formula is identical, it is just fed the state of a different
% plant, so the commanded values come out different. This is what the real
% robot would actually be asked to do.
fprintf('(running the exact-kinematics model for comparison ...)\n');
x0_rad = [P.x0(1); P.x0(2); deg2rad(P.x0(3)); deg2rad(P.x0(4))];
[t_ex, ~, z_ex, u1_ex, u2_ex] = simulate_case_physical( ...
    x0_rad, P.d0, P.l1, P.alpha, P.gamma, P.epsilon, P.omega, ...
    P.k12, P.k112, P.t_start, P.t_end, P.dt, P.h_fun);

% ---------------- state trajectories ----------------------------------------
figure('Name', 'State trajectories (1 trailer)', 'NumberTitle', 'off');
labels = {'z_1', 'z_2', 'z_3', 'z_4'};
for i = 1:4
    subplot(4,1,i)
    plot(t_sim, z_sim(:,i), 'Color', col, 'LineWidth', 1.2);
    grid on; ylabel(labels{i});
    if i == 1
        title('State trajectories');
    end
    if i == 4
        xlabel('t [s]');
    end
end

% ---------------- control inputs ---------------------------------------------
figure('Name', 'Control inputs (1 trailer)', 'NumberTitle', 'off');
subplot(2,1,1)
plot(t_sim, u1_sim, 'Color', col,    'LineWidth', 1.0); hold on
plot(t_ex,  u1_ex,  'Color', col_ex, 'LineWidth', 1.0);
yline(U1_MAX, 'r--'); yline(-U1_MAX, 'r--');
grid on; ylabel('u_1 [m/s]');
title('Control inputs');

subplot(2,1,2)
plot(t_sim, u2_sim, 'Color', col,    'LineWidth', 1.0); hold on
plot(t_ex,  u2_ex,  'Color', col_ex, 'LineWidth', 1.0);
yline(U2_MAX, 'r--'); yline(-U2_MAX, 'r--');
grid on; xlabel('t [s]'); ylabel('u_2 [rad/s]');

% ---------------- convergence -------------------------------------------------
figure('Name', 'Convergence (1 trailer)', 'NumberTitle', 'off');
semilogy(t_sim, vecnorm(z_sim, 2, 2), 'Color', col, 'LineWidth', 1.5);
grid on; xlabel('t [s]'); ylabel('||z(t)||');
title('Convergence');

legend('approximation', 'exact kinematics', 'limit', 'Location', 'best');

r1 = max(abs(u1_sim))/U1_MAX;
r2 = max(abs(u2_sim))/U2_MAX;
r1e = max(abs(u1_ex))/U1_MAX;
r2e = max(abs(u2_ex))/U2_MAX;
fprintf('\n  max |u1| = %10.2f m/s    limit %.2f    -> %8.1f x over\n', max(abs(u1_sim)), U1_MAX, r1);
fprintf('  max |u2| = %10.2f rad/s  limit %.2f    -> %8.1f x over\n', max(abs(u2_sim)), U2_MAX, r2);
fprintf('\n  --- same controller, EXACT kinematics (what the robot does) ---\n');
fprintf('  max |u1| = %10.2f m/s    limit %.2f    -> %8.1f x over\n', max(abs(u1_ex)), U1_MAX, r1e);
fprintf('  max |u2| = %10.2f rad/s  limit %.2f    -> %8.1f x over\n', max(abs(u2_ex)), U2_MAX, r2e);
fprintf('  final ||z|| : approximation %.3e,  exact %.3e\n', ...
        norm(z_sim(end,:)), norm(z_ex(end,:)));
if r1e > 1 || r2e > 1
    fprintf(2, '  *** the EXACT model exceeds a limit -- see OBALKA_1PRIVES.md ***\n');
end
