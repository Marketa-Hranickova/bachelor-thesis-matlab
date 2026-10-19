% CONVERGENCE_MAP  Does (epsilon, gamma) converge, overshoot, or diverge?
% ONE-TRAILER (r=3, 4-state) version -- see "2 trailers/convergence_map.m"
% for the full explanation of the method and the linear-bound reference
% lines. Same idea, just the 4-state system and 2 multipliers (k12,k112).
% Uses the Lyapunov function from params.m, so the map is for the
% controller actually in use.
%
% Run:  convergence_map

clc; close all;

P = params();
alpha = P.alpha;
k12 = P.k12;  k112 = P.k112;
z0 = P.z0;  t_start = P.t_start;

h_fun  = P.h_fun;    % Lyapunov function lives in params.m

t_end_map = 8;
dt_map    = 1e-3;

eps_list   = logspace(log10(0.05), log10(16), 12);
gamma_list = logspace(log10(0.01), log10(3),  12);

n0 = norm(z0);
verdict = zeros(numel(eps_list), numel(gamma_list));   % 1=converges 2=marginal 3=diverges

fprintf('Sweeping %d x %d = %d combinations, t_end = %g s ...\n', ...
        numel(eps_list), numel(gamma_list), numel(eps_list)*numel(gamma_list), t_end_map);
tic
for i = 1:numel(eps_list)
    Q       = set_epsilon(P, eps_list(i));   % epsilon AND omega together
    eps_i   = Q.epsilon;
    omega_i = Q.omega;
    for j = 1:numel(gamma_list)
        g_j = gamma_list(j);
        try
            [~, z_sim] = simulate_case(z0, alpha, g_j, ...
                eps_i, omega_i, k12, k112, ...
                t_start, t_end_map, dt_map, h_fun);
            zf = z_sim(end, :);
        catch
            zf = [NaN NaN NaN NaN];
        end

        if any(~isfinite(zf))
            verdict(i, j) = 3;
        else
            ratio = norm(zf) / n0;
            if ratio < 0.1
                verdict(i, j) = 1;
            elseif ratio < 2
                verdict(i, j) = 2;
            else
                verdict(i, j) = 3;
            end
        end
    end
end
fprintf('done in %.1f s\n', toc);

[E, G] = meshgrid(eps_list, gamma_list);
E = E.';  G = G.';

colors = [0.20 0.65 0.30
          0.90 0.70 0.10
          0.80 0.20 0.20];

figure('Name', 'Convergence map (epsilon, gamma) -- 1 trailer', 'NumberTitle', 'off');
hold on

h = gobjects(1, 3);
labels = {'converges (||z_{end}|| < 10% ||z_0||)', 'marginal', 'diverges / blows up'};
for k = 1:3
    mask = verdict == k;
    if any(mask(:))
        h(k) = scatter(E(mask), G(mask), 70, colors(k,:), 'filled', ...
                       'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    else
        h(k) = scatter(nan, nan, 70, colors(k,:), 'filled');
    end
end

eps_ref = logspace(log10(min(eps_list)), log10(max(eps_list)), 100);
r1 = plot(eps_ref, 1./eps_ref, 'k--', 'LineWidth', 1.2);
r2 = plot(eps_ref, 2./eps_ref, 'k:',  'LineWidth', 1.0);

set(gca, 'XScale', 'log', 'YScale', 'log');
grid on
xlabel('\epsilon  [s]'); ylabel('\gamma');
legend([h, r1, r2], [labels, ...
       {'\epsilon\gamma = 1  (proof bound: exact hit)', ...
        '\epsilon\gamma = 2  (proof bound: predicted breakdown)'}], ...
       'Location', 'eastoutside');
title(sprintf('True convergence vs. the linear bound''s prediction   (1 trailer, k = %d, %d)', ...
              k12, k112));

fprintf('\n=== Done. ===\n');
