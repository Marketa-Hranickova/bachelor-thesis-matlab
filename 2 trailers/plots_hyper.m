function plots_hyper(opts)
%PLOTS_HYPER  Plots for the folder "hyperexponential 2 trailers" (task from Victoria, 23 Oct 2026),
%   in the style of plots_o.m (folder original_system).
%
%   EXACT kinematics of the car with two trailers, parameters of the paper
%   (Zuyev & Grushkovskaya), geometry (26): d0 = d1 = 0.1, l1 = l2 = 1,
%   eps = 0.1 s, gamma = 1, k = (8,6,5), NO velocity limits.
%   Start of the paper: z0 = (1,-2,2,-1,1), x0 from z_to_x.
%   Blue solid: h = 1.   Green dashed: h = |log(||z||^2)|  (same colours as plots_o.m).
%
%   Figures (shown, and saved as PNG into ../hyperexponential 2 trailers/):
%     x_trajectories   x_C0, y_C0, phi, theta_1, theta_2 against t
%     control_inputs   u_1(t), u_2(t)
%     norm_z           ||z(t)|| / ||z_0||
%
%   plots_hyper                         full run (16 s for both h)
%   plots_hyper(struct('Tscale',0.2))   shorter runs (quick check)
%   plots_hyper(struct('show',false))   do not display the figures, only save them

    if nargin < 1, opts = struct(); end
    if ~isfield(opts, 'Tscale'), opts.Tscale = 1; end
    if ~isfield(opts, 'show'),   opts.show = true; end
    here = fileparts(mfilename('fullpath'));
    if ~isfield(opts, 'outdir'), opts.outdir = fullfile(here, '..', 'hyperexponential 2 trailers'); end
    if ~exist(opts.outdir, 'dir'), mkdir(opts.outdir); end
    vis = 'on';  if ~opts.show, vis = 'off'; end
    setappdata(0, 'ctrl_phase', [0 0 0]);
    Tr = @(name, varargin) trailer2(name, varargin{:});

    g = struct('d0', 0.1, 'd1', 0.1, 'l1', 1, 'l2', 1);  c = [0.10495 0.0024 0.0112];
    ep = 0.1;  om = 2*pi/ep;  gam = 1;  k = [8 6 5];
    dt = min(0.02, 2*pi/(max(k)*om)/40);
    z0 = [1; -2; 2; -1; 1];  r0 = norm(z0);
    x0 = Tr('z_to_x', z0, g.d0, g.d1, g.l1, g.l2);
    hf = {@(z) 1, @(z) abs(log(sum(z.^2)))};
    T  = [16 16] * opts.Tscale;

    for hi = 1:2
        fprintf('run %d of 2 (T = %g s) ...\n', hi, T(hi));
        [t, x, z, u1, u2] = Tr('simulate_case_physical', x0, g.d0, g.d1, g.l1, g.l2, c(1), c(2), c(3), ...
                               gam, ep, om, k(1), k(2), k(3), 0, T(hi), dt, hf{hi});
        D(hi).t = t;  D(hi).x = x;  D(hi).nz = sqrt(sum(z.^2, 2))/r0;  D(hi).u1 = u1;  D(hi).u2 = u2;
    end

    % ---- plot style (as in plots_o.m) ------------------------------------------
    c1_col = [0.00 0.20 0.70];
    c2_col = [0.00 0.70 0.40];
    x_labels = {'x_{C_0}', 'y_{C_0}', '\phi', '\theta_1', '\theta_2'};
    leg = {'h = 1', 'h = |log(||z||^2)|'};

    % ---- x trajectories -----------------------------------------------------------
    f = figure('Name', 'x Trajectories', 'NumberTitle', 'off', 'Visible', vis);
    for kk = 1:5
        subplot(5,1,kk)
        plot(D(1).t, D(1).x(:,kk), 'Color', c1_col, 'LineWidth', 1.5, 'LineStyle', '-'); hold on;
        plot(D(2).t, D(2).x(:,kk), 'Color', c2_col, 'LineWidth', 1.5, 'LineStyle', '--');
        grid on;  ylabel(x_labels{kk});
        if kk == 1
            title('x-Trajectories (physical state) - paper parameters (26), start z_0 = (1,-2,2,-1,1)');
            legend(leg, 'Location', 'best');
        end
        if kk == 5, xlabel('Time [s]'); end
    end
    save_png(f, fullfile(opts.outdir, 'x_trajectories.png'));

    % ---- control inputs ---------------------------------------------------------------
    f = figure('Name', 'Control Inputs', 'NumberTitle', 'off', 'Visible', vis);
    subplot(2,1,1)
    plot(D(1).t, D(1).u1, 'Color', c1_col, 'LineWidth', 1.2, 'LineStyle', '-'); hold on;
    plot(D(2).t, D(2).u1, 'Color', c2_col, 'LineWidth', 1.2, 'LineStyle', '--');
    grid on;  ylabel('u_1(t)');
    title('Control Inputs - paper parameters (26), start z_0 = (1,-2,2,-1,1)');
    legend(leg, 'Location', 'best');
    subplot(2,1,2)
    plot(D(1).t, D(1).u2, 'Color', c1_col, 'LineWidth', 1.2, 'LineStyle', '-'); hold on;
    plot(D(2).t, D(2).u2, 'Color', c2_col, 'LineWidth', 1.2, 'LineStyle', '--');
    grid on;  ylabel('u_2(t)');  xlabel('Time [s]');
    legend(leg, 'Location', 'best');
    save_png(f, fullfile(opts.outdir, 'control_inputs.png'));

    % ---- norm of z --------------------------------------------------------------------------
    f = figure('Name', 'Norm of z', 'NumberTitle', 'off', 'Visible', vis);
    plot(D(1).t, D(1).nz, 'Color', c1_col, 'LineWidth', 1.5, 'LineStyle', '-'); hold on;
    plot(D(2).t, D(2).nz, 'Color', c2_col, 'LineWidth', 1.5, 'LineStyle', '--');
    grid on;  xlabel('Time [s]');  ylabel('||z(t)|| / ||z_0||');
    title('Norm of derived z(t) - paper parameters (26), start z_0 = (1,-2,2,-1,1)');
    legend(leg, 'Location', 'best');
    save_png(f, fullfile(opts.outdir, 'norm_z.png'));

    fprintf('PNG files written to %s\n', opts.outdir);
end

function save_png(f, file)
% Save figure f as PNG without the axes toolbar (MATLAB R2018b+ warns otherwise).
    for ax = findall(f, 'Type', 'axes')'
        try, ax.Toolbar.Visible = 'off'; catch, end      % older MATLAB / Octave: no toolbar, nothing to do
    end
    print(f, file, '-dpng', '-r150');
end
