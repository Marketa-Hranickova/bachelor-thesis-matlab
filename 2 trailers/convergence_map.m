

clc; close all;

P = params();
alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;
k12 = P.k12;  k112 = P.k112;  k1112 = P.k1112;
z0 = P.z0;  t_start = P.t_start;

h_fun  = @(z) 1;
g1_fun = @(z, alpha, kappa, theta) [1; 0; -z(2); -z(3); kappa*z(2)^3 - z(4)];
g2_fun = @(z, alpha, kappa, theta) [0; 1; 0; -alpha*z(2)^2; -theta*z(2)*z(3)];

% ---------------- sweep settings (fast, short horizon) --------------------
t_end_map = 8;          % seconds -- enough periods to see the trend
dt_map    = 1e-3;       % coarser than calculations.m's 1e-4, for speed

eps_list   = logspace(log10(0.05), log10(16), 12);
gamma_list = logspace(log10(0.01), log10(3),  12);

n0 = norm(z0);
% verdict: 1 = converges, 2 = marginal, 3 = diverges / blows up
verdict = zeros(numel(eps_list), numel(gamma_list));

fprintf('Sweeping %d x %d = %d combinations, t_end = %g s ...\n', ...
        numel(eps_list), numel(gamma_list), numel(eps_list)*numel(gamma_list), t_end_map);
tic
for i = 1:numel(eps_list)
    eps_i   = eps_list(i);
    omega_i = 2*pi/eps_i;
    for j = 1:numel(gamma_list)
        g_j = gamma_list(j);
        try
            [~, z_sim] = simulate_case(z0, alpha, kappa, theta, g_j, ...
                eps_i, omega_i, k12, k112, k1112, ...
                t_start, t_end_map, dt_map, h_fun, g1_fun, g2_fun);
            zf = z_sim(end, :);
        catch
            zf = [NaN NaN NaN NaN NaN];   % solver choked -> treat as divergence
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

% ---------------- numeric check: real boundary vs. eps*gamma = 1/2 ---------
fprintf('\n%-8s %14s %10s   %14s %10s\n', ...
        'gamma', 'last conv eps', 'eps*gamma', 'first div eps', 'eps*gamma');
fprintf('%s\n', repmat('-', 1, 62));
for j = 1:numel(gamma_list)
    g_j = gamma_list(j);
    conv_idx = find(verdict(:, j) == 1);
    div_idx  = find(verdict(:, j) == 3);

    if isempty(conv_idx)
        last_conv_str = '     --   ';
        last_conv_eg  = NaN;
    else
        last_conv_eps = eps_list(max(conv_idx));
        last_conv_str = sprintf('%10.4f', last_conv_eps);
        last_conv_eg  = last_conv_eps * g_j;
    end

    if isempty(div_idx)
        first_div_str = '     --   ';
        first_div_eg  = NaN;
    else
        first_div_eps = eps_list(min(div_idx));
        first_div_str = sprintf('%10.4f', first_div_eps);
        first_div_eg  = first_div_eps * g_j;
    end

    fprintf('%-8.4f %14s %10.3f   %14s %10.3f\n', ...
            g_j, last_conv_str, last_conv_eg, first_div_str, first_div_eg);
end
fprintf('\n(compare the eps*gamma columns to the proof bound''s 1 and 2)\n');

% ---------------- plot ------------------------------------------------------
[E, G] = meshgrid(eps_list, gamma_list);
E = E.';  G = G.';

colors = [0.20 0.65 0.30      % converges -- green
          0.90 0.70 0.10      % marginal  -- amber
          0.80 0.20 0.20];    % diverges  -- red

figure('Name', 'Convergence map (epsilon, gamma)', 'NumberTitle', 'off');
hold on

h = gobjects(1, 3);
labels = {'converges (||z_{end}|| < 10% ||z_0||)', 'marginal', 'diverges / blows up'};
for k = 1:3
    mask = verdict == k;
    if any(mask(:))
        h(k) = scatter(E(mask), G(mask), 70, colors(k,:), 'filled', ...
                       'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    else
        h(k) = scatter(nan, nan, 70, colors(k,:), 'filled');  % keep legend entry
    end
end

% reference curves from the LINEARISED proof bound
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
title(sprintf('True convergence vs. the linear bound''s prediction   (k = %d, %d, %d)', ...
              k12, k112, k1112));

fprintf('\nColored points = true nonlinear simulation result.\n');
fprintf('Dashed/dotted lines = where the LINEAR proof bound predicts the boundary.\n');
fprintf('If the real boundary sits noticeably below eps*gamma = 1, the nonlinear\n');
fprintf('terms are biting earlier than the linearised proof suggests -- worth a note.\n');
fprintf('\n=== Done. ===\n');
