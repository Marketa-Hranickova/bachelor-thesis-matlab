function thesis_figures_compare(fig_dir)
%THESIS_FIGURES_COMPARE  r = 2, 3, 4 side by side (thesis Chapter 6, "degree
%   of nonholonomy"). Run AFTER thesis_figures in "1 trailer" and in
%   "2 trailers": it reads the feas_1tr.mat and feas_2tr.mat they wrote.
%
%   thesis_figures_compare            reads 1 trailer/figures and 2 trailers/figures
%   thesis_figures_compare(folder)    reads both .mat files from one folder
%
%   r = 3 and r = 4 are taken from the .mat files, i.e. from the same code
%   as the rest of the study. r = 2 (car alone, no trailer: only the bracket
%   [f1,f2] is needed) has no code of its own in this study, so it is
%   computed here with the same recipe, restricted to its terms:
%
%       u1 = a1 + T12 cos(k12 w t),    u2 = a2 + T12 sign(a12) sin(k12 w t)
%       a1 = -gamma z1,  a2 = -gamma z2,  a12 = -gamma (z1 z2 + z3)
%       T12 = sqrt(2 k12 w |a12|),   k12 = 1
%
%   from the first three privileged coordinates of the same start pose.
%   Quadratic V (h = 1), as in the sweeps of thesis_figures.

    here = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(fig_dir)
        dirA = fullfile(here, '1 trailer', 'figures');   % where thesis_figures of each folder writes by default
        dirB = fullfile(here, '2 trailers', 'figures');
    else
        dirA = fig_dir;  dirB = fig_dir;
    end
    A = load(fullfile(dirA, 'feas_1tr.mat'));
    B = load(fullfile(dirB, 'feas_2tr.mat'));
    assert(isequal(A.eps_list, B.eps_list), 'the two runs used different epsilon lists');
    eps_list = A.eps_list;  U1 = A.u1max;  U2 = A.u2max;  cap = A.cap;

    % ---------------- r = 2 ----------------------------------------------------
    z = A.z0(1:3);  k12 = 1;
    g2 = zeros(size(eps_list));
    for i = 1:numel(eps_list)
        w  = 2*pi/eps_list(i);
        lo = 1e-16;  hi = 10;
        for it = 1:200
            gam = sqrt(lo*hi);
            a1 = -gam*z(1);  a2 = -gam*z(2);  a12 = -gam*(z(1)*z(2) + z(3));
            T12 = sqrt(2*k12*w*abs(a12));
            if abs(a1) + T12 <= U1 && abs(a2) + T12 <= U2
                lo = gam;
            else
                hi = gam;
            end
        end
        g2(i) = lo;
    end
    g2_eff = min(g2, cap./eps_list);
    t2 = log(10) ./ g2_eff;

    % ---------------- table for the text -------------------------------------
    fprintf('\n   eps [s]   t(10%%) r=2 [min]   r=3 [min]   r=4 [h]   r=4 / r=3\n');
    for i = 1:numel(eps_list)
        fprintf('  %6.2f   %14.3f   %10.2f   %8.2f   %9.0f\n', eps_list(i), ...
                t2(i)/60, A.t10(i)/60, B.t10(i)/3600, B.t10(i)/A.t10(i));
    end
    fprintf('\n   eps [s]   gamma_amp  r=2        r=3        r=4\n');
    for i = 1:numel(eps_list)
        fprintf('  %6.2f   %10.3e %10.3e %10.3e\n', eps_list(i), g2(i), A.g_amp(i), B.g_amp(i));
    end
    % local slopes of gamma_amp(epsilon) at small epsilon (log-log)
    sl = @(g) diff(log(g(1:3))) ./ diff(log(eps_list(1:3)));
    fprintf('\nlog-log slope of gamma_amp over eps = %g..%g s:  r=2 %.2f %.2f   r=3 %.2f %.2f   r=4 %.2f %.2f\n', ...
            eps_list(1), eps_list(3), sl(g2), sl(A.g_amp), sl(B.g_amp));

    % ---------------- figure --------------------------------------------------
    c2 = [0.47 0.67 0.19];  c3 = [0.00 0.45 0.74];  c4 = [0.85 0.33 0.10];
    f = figure('Color', 'w');
    set(f, 'DefaultAxesFontSize', 9, 'DefaultTextFontSize', 9);

    subplot(1,2,1)
    loglog(eps_list, g2,        '^-', 'Color', c2, 'LineWidth', 1.5); hold on
    loglog(eps_list, A.g_amp,   'o-', 'Color', c3, 'LineWidth', 1.5);
    loglog(eps_list, B.g_amp,   's-', 'Color', c4, 'LineWidth', 1.5);
    loglog(eps_list, cap./eps_list, 'k:', 'LineWidth', 1.0);
    grid on; xlabel('\epsilon  [s]'); ylabel('\gamma_{amp}');

    subplot(1,2,2)
    loglog(eps_list, t2/3600,      '^-', 'Color', c2, 'LineWidth', 1.5); hold on
    loglog(eps_list, A.t10/3600,   'o-', 'Color', c3, 'LineWidth', 1.5);
    loglog(eps_list, B.t10/3600,   's-', 'Color', c4, 'LineWidth', 1.5);
    yl = ylim;
    patch([4 max(eps_list) max(eps_list) 4], [yl(1) yl(1) yl(2) yl(2)], c4, ...
          'FaceAlpha', 0.08, 'EdgeColor', 'none');
    ylim(yl);
    grid on; xlabel('\epsilon  [s]'); ylabel('time to 10 % of ||z_0||  [h]');
    legend('r = 2', 'r = 3', 'r = 4', 'Location', 'northeast', 'Box', 'off');

    set(f, 'Units', 'centimeters', 'Position', [2 2 15 6.5], ...
           'PaperUnits', 'centimeters', 'PaperSize', [15 6.5], 'PaperPosition', [0 0 15 6.5]);
    print(f, fullfile(fig_dir, 'feas_cmp_eps'), '-dpdf', '-painters');
    close(f);
    fprintf('\nfeas_cmp_eps.pdf written to %s\n', fig_dir);
end
