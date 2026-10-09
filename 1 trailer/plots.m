% plots.m -- figures for the ONE-TRAILER (r=3, 4-state) simulation.
% Run calculations.m first (it leaves t_sim, z_sim, u1_sim, u2_sim, P in
% the workspace).
%
% Run:  plots

P = params();
U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

col = [0.00 0.20 0.70];

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
plot(t_sim, u1_sim, 'Color', col, 'LineWidth', 1.0); hold on
yline(U1_MAX, 'r--'); yline(-U1_MAX, 'r--');
grid on; ylabel('u_1 [m/s]');
title('Control inputs');

subplot(2,1,2)
plot(t_sim, u2_sim, 'Color', col, 'LineWidth', 1.0); hold on
yline(U2_MAX, 'r--'); yline(-U2_MAX, 'r--');
grid on; xlabel('t [s]'); ylabel('u_2 [rad/s]');

% ---------------- convergence -------------------------------------------------
figure('Name', 'Convergence (1 trailer)', 'NumberTitle', 'off');
semilogy(t_sim, vecnorm(z_sim, 2, 2), 'Color', col, 'LineWidth', 1.5);
grid on; xlabel('t [s]'); ylabel('||z(t)||');
title('Convergence');

r1 = max(abs(u1_sim))/U1_MAX;
r2 = max(abs(u2_sim))/U2_MAX;
fprintf('\n  max |u1| = %10.2f m/s    limit %.2f    -> %8.1f x over\n', max(abs(u1_sim)), U1_MAX, r1);
fprintf('  max |u2| = %10.2f rad/s  limit %.2f    -> %8.1f x over\n', max(abs(u2_sim)), U2_MAX, r2);
