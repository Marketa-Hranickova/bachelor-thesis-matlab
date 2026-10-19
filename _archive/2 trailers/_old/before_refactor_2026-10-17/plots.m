% plots.m -- figures for the single simulated case.
%
% Run this; it calls calculations.m, which calls params.m.
% Three figures: states, controls, convergence.

run('calculations.m');

% -------------------------------------------------------------- style
col_main  = [0.00 0.20 0.70];    % blue
col_alt   = [0.85 0.33 0.10];    % orange
col_limit = [0.45 0.45 0.45];    % grey, for the hardware limits

U1_MAX = 0.22;                   % TurtleBot3 Burger  [m/s]
U2_MAX = 2.84;                   % TurtleBot3 Burger  [rad/s]

state_labels = {'z_1', 'z_2', 'z_3', 'z_4', 'z_5'};
geom_str = sprintf('d_0=%.3f, l_1=%.3f, d_1=%.3f, l_2=%.3f m', d0, l1, d1, l2);

% =========================================================================
%  FIGURE 1 -- state trajectories
% =========================================================================
figure('Name', 'State trajectories', 'NumberTitle', 'off');
for k = 1:5
    subplot(5, 1, k)
    plot(t_sim, z_sim(:, k), 'Color', col_main, 'LineWidth', 1.4);
    grid on
    ylabel(state_labels{k});
    if k == 1
        title(sprintf('State trajectories  (%s)', geom_str));
    end
    if k == 5
        xlabel('Time [s]');
    end
end

% =========================================================================
%  FIGURE 2 -- control inputs, against the hardware limits
% =========================================================================
figure('Name', 'Control inputs', 'NumberTitle', 'off');

subplot(2, 1, 1)
plot(t_sim, u1_sim, 'Color', col_main, 'LineWidth', 1.0); hold on
yline( U1_MAX, '--', 'Color', col_limit, 'LineWidth', 1.0);
yline(-U1_MAX, '--', 'Color', col_limit, 'LineWidth', 1.0);
grid on
ylabel('u_1(t)  [m/s]');
title(sprintf('Control inputs   (dashed: TurtleBot3 limits %.2f m/s, %.2f rad/s)', ...
              U1_MAX, U2_MAX));

subplot(2, 1, 2)
plot(t_sim, u2_sim, 'Color', col_alt, 'LineWidth', 1.0); hold on
yline( U2_MAX, '--', 'Color', col_limit, 'LineWidth', 1.0);
yline(-U2_MAX, '--', 'Color', col_limit, 'LineWidth', 1.0);
grid on
ylabel('u_2(t)  [rad/s]');
xlabel('Time [s]');

% =========================================================================
%  FIGURE 3 -- convergence
% =========================================================================
% ||z(t)|| on a log axis: exponential convergence shows as a straight line.
figure('Name', 'Convergence', 'NumberTitle', 'off');
semilogy(t_sim, max(sqrt(V_sim), eps), 'Color', col_main, 'LineWidth', 1.5);
grid on
xlabel('Time [s]');
ylabel('||z(t)||');
title('Convergence of the state norm');

% =========================================================================
%  How far outside the robot's limits is this?
% =========================================================================
r1 = max(abs(u1_sim))/U1_MAX;
r2 = max(abs(u2_sim))/U2_MAX;
fprintf('\n========================================\n')
fprintf('  Control demand vs TurtleBot3 Burger\n')
fprintf('========================================\n')
fprintf('  max |u1| = %10.2f m/s    limit %.2f    -> %8.0f x over\n', ...
        max(abs(u1_sim)), U1_MAX, r1);
fprintf('  max |u2| = %10.2f rad/s  limit %.2f    -> %8.0f x over\n', ...
        max(abs(u2_sim)), U2_MAX, r2);
if r1 > 1 || r2 > 1
    fprintf('  -> the commanded controls do NOT fit on this robot.\n');
else
    fprintf('  -> the commanded controls fit within the limits.\n');
end
