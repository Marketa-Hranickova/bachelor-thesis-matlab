% PLOTS  One closed-loop run of the approximate (z) model, with figures.
%
%   Run this to see a single simulation: the states, the controls against
%   the TurtleBot limits, and the convergence of ||z||.
%   Everything comes from params.m -- change numbers there.
%
%   This is the nilpotent APPROXIMATION, the system the controller is
%   designed on. What the robot really does is in compare_physical.m.

P = params();
print_setup(P);

fprintf('\nsimulating (approximate model) ...\n');
tic
out = simulate(P);
fprintf('  done in %.1f s, %d samples\n', toc, numel(out.t));

t   = out.t;
z   = out.z;
nz  = sqrt(sum(z.^2, 2));                  % ||z(t)||

% ---------------------------------------------------------------- style
col_main  = [0.00 0.20 0.70];              % blue
col_alt   = [0.85 0.33 0.10];              % orange
col_limit = [0.45 0.45 0.45];              % grey, for the hardware limits
U1 = P.u1_max;
U2 = P.u2_max;
geom_str = sprintf('d_0=%.3f, l_1=%.3f, d_1=%.3f, l_2=%.3f m', P.d0, P.l1, P.d1, P.l2);

% =========================================================================
%  FIGURE 1 -- state trajectories
% =========================================================================
figure('Name', 'State trajectories', 'NumberTitle', 'off');
for k = 1:5
    subplot(5, 1, k)
    plot(t, z(:, k), 'Color', col_main, 'LineWidth', 1.4);
    grid on
    ylabel(sprintf('z_%d', k));
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
plot(t, out.u1, 'Color', col_main, 'LineWidth', 1.0); hold on
yline( U1, '--', 'Color', col_limit, 'LineWidth', 1.0);
yline(-U1, '--', 'Color', col_limit, 'LineWidth', 1.0);
grid on
ylabel('u_1(t)  [m/s]');
title(sprintf('Control inputs   (dashed: TurtleBot3 limits %.2f m/s, %.2f rad/s)', U1, U2));

subplot(2, 1, 2)
plot(t, out.u2, 'Color', col_alt, 'LineWidth', 1.0); hold on
yline( U2, '--', 'Color', col_limit, 'LineWidth', 1.0);
yline(-U2, '--', 'Color', col_limit, 'LineWidth', 1.0);
grid on
ylabel('u_2(t)  [rad/s]');
xlabel('Time [s]');

% =========================================================================
%  FIGURE 3 -- convergence
% =========================================================================
% ||z(t)|| on a log axis: exponential convergence shows as a straight line.
figure('Name', 'Convergence', 'NumberTitle', 'off');
semilogy(t, max(nz, eps), 'Color', col_main, 'LineWidth', 1.5);
grid on
xlabel('Time [s]');
ylabel('||z(t)||');
title(sprintf('Convergence of the state norm   (V = %s)', P.lyap_name));

% =========================================================================
%  Summary: result and how close the controls come to the robot's limits
% =========================================================================
m1 = max(abs(out.u1));
m2 = max(abs(out.u2));
fprintf('\n========================================\n')
fprintf('  Result\n')
fprintf('========================================\n')
fprintf('  ||z||   : %.6f  ->  %.6f   (%.2f %% of start)\n', ...
        nz(1), nz(end), 100*nz(end)/nz(1));
fprintf('  final z : [%.6f %.6f %.6f %.6f %.6f]\n', z(end, :));

fprintf('\n========================================\n')
fprintf('  Control demand vs TurtleBot3 Burger\n')
fprintf('========================================\n')
fprintf('  max |u1| = %8.4f m/s    limit %.2f    -> %.2f x the limit\n', m1, U1, m1/U1);
fprintf('  max |u2| = %8.4f rad/s  limit %.2f    -> %.2f x the limit\n', m2, U2, m2/U2);
if m1 > U1 || m2 > U2
    fprintf('  -> the commanded controls do NOT fit on this robot.\n');
else
    fprintf('  -> the commanded controls fit within the limits.\n');
end
fprintf('  (approximate model only -- run compare_physical for the exact one)\n');
