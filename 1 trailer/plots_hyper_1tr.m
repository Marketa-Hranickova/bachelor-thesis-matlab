function plots_hyper_1tr(opts)
%PLOTS_HYPER_1TR  ONE trailer: h = |log(||z||^2)| (hyperexponential try) against Victoria's log V (c = 10),
%   in the style of plots_o.m (folder original_system) and of ../2 trailers/plots_hyper.m.
%
%   EXACT kinematics, robot geometry (params.m: d0 = 0.115, l1 = 0.082), eps = 5 s, k = (3,2),
%   start x0 = (-0.2, -0.2, 45 deg, 45 deg). Both runs within the velocity limits:
%     blue solid   : Victoria's log V, h = 1/(1 + c ||z||^2), c = 10, gamma from params.m
%     green dashed : h = |log(||z||^2)|, gamma = largest admissible along the ray (gamma_ray, cap 0.8)
%
%   Figures (shown, and saved as PNG into ../hyperexponential 1 trailer/):
%     x_trajectories   x_C0, y_C0, phi, theta_1 against t
%     V_functions      V(z) = ||z||^2  and  V_2(z) = ||z||^2 |log(||z||^2)|   (as in plots_o.m)
%     control_inputs   u_1(t), u_2(t)
%     norm_z           ||z(t)|| / ||z_0||
%
%   plots_hyper_1tr                        full run (500 s)
%   plots_hyper_1tr(struct('Tscale',0.2))  shorter runs (quick check)
%   plots_hyper_1tr(struct('show',false))  do not display, only save
%
%   Call it from the folder "1 trailer" (not together with "2 trailers" on the path).

    if nargin < 1, opts = struct(); end
    if ~isfield(opts, 'Tscale'), opts.Tscale = 1; end
    if ~isfield(opts, 'show'),   opts.show = true; end
    here = fileparts(mfilename('fullpath'));
    if ~isfield(opts, 'outdir'), opts.outdir = fullfile(here, '..', 'hyperexponential 1 trailer'); end
    if ~exist(opts.outdir, 'dir'), mkdir(opts.outdir); end
    vis = 'on';  if ~opts.show, vis = 'off'; end
    Tr = @(name, varargin) trailer1(name, varargin{:});

    P  = params();
    x0 = [P.x0(1:2); deg2rad(P.x0(3:4))];                 % radians for the simulator
    z0 = Tr('x_to_z', x0, P.d0, P.l1);  r0 = norm(z0);
    h_log = @(z) abs(log(sum(z.^2)));
    gam_log = Tr('gamma_ray', z0, P.alpha, P.omega, P.epsilon, P.k12, P.k112, ...
                 P.u1_max, P.u2_max, h_log, P.cap_product);
    hf  = {P.h_fun, h_log};
    gam = [P.gamma, gam_log];
    fprintf('gamma: log V (c = 10) %.4e, |log| %.4e\n', gam);
    T  = 500 * opts.Tscale;
    dt = min(0.02, 2*pi/(max(P.k12, P.k112)*P.omega)/40);

    for hi = 1:2
        fprintf('run %d of 2 (T = %g s) ...\n', hi, T);
        [t, x, z, u1, u2] = Tr('simulate_case_physical', x0, P.d0, P.l1, P.alpha, gam(hi), P.epsilon, P.omega, ...
                               P.k12, P.k112, 0, T, dt, hf{hi});
        D(hi).t = t;  D(hi).x = x;  D(hi).z = z;  D(hi).u1 = u1;  D(hi).u2 = u2;
        D(hi).V  = sum(z.^2, 2);
        D(hi).V2 = D(hi).V .* abs(log(D(hi).V));  D(hi).V2(D(hi).V == 0) = 0;
        D(hi).nz = sqrt(D(hi).V)/r0;
    end

    % ---- plot style (as in plots_o.m) ------------------------------------------
    c1_col = [0.00 0.20 0.70];
    c2_col = [0.00 0.70 0.40];
    x_labels = {'x_{C_0}', 'y_{C_0}', '\phi', '\theta_1'};
    leg = {'log V, c = 10', 'h = |log(||z||^2)|'};
    ttl = 'one trailer, robot geometry, \epsilon = 5 s, k = (3,2)';

    % ---- x trajectories -----------------------------------------------------------
    f = figure('Name', 'x Trajectories', 'NumberTitle', 'off', 'Visible', vis);
    for kk = 1:4
        subplot(4,1,kk)
        plot(D(1).t, D(1).x(:,kk), 'Color', c1_col, 'LineWidth', 1.5, 'LineStyle', '-'); hold on;
        plot(D(2).t, D(2).x(:,kk), 'Color', c2_col, 'LineWidth', 1.5, 'LineStyle', '--');
        grid on;  ylabel(x_labels{kk});
        if kk == 1, title(['x-Trajectories (physical state) - ' ttl]);  legend(leg, 'Location', 'best'); end
        if kk == 4, xlabel('Time [s]'); end
    end
    save_png(f, fullfile(opts.outdir, 'x_trajectories.png'));

    % ---- V functions ------------------------------------------------------------------
    f = figure('Name', 'V(z) and V2(z)', 'NumberTitle', 'off', 'Visible', vis);
    subplot(2,1,1)
    plot(D(1).t, D(1).V, 'Color', c1_col, 'LineWidth', 1.5, 'LineStyle', '-'); hold on;
    plot(D(2).t, D(2).V, 'Color', c2_col, 'LineWidth', 1.5, 'LineStyle', '--');
    grid on;  ylabel('V(z) = ||z||^2');  title(['V(z) and V_2(z) - ' ttl]);  legend(leg, 'Location', 'best');
    subplot(2,1,2)
    plot(D(1).t, D(1).V2, 'Color', c1_col, 'LineWidth', 1.5, 'LineStyle', '-'); hold on;
    plot(D(2).t, D(2).V2, 'Color', c2_col, 'LineWidth', 1.5, 'LineStyle', '--');
    grid on;  ylabel('V_2(z) = ||z||^2 * |log(||z||^2)|');  xlabel('Time [s]');  legend(leg, 'Location', 'best');
    save_png(f, fullfile(opts.outdir, 'V_functions.png'));

    % ---- control inputs ---------------------------------------------------------------
    f = figure('Name', 'Control Inputs', 'NumberTitle', 'off', 'Visible', vis);
    subplot(2,1,1)
    plot(D(1).t, D(1).u1, 'Color', c1_col, 'LineWidth', 1.2, 'LineStyle', '-'); hold on;
    plot(D(2).t, D(2).u1, 'Color', c2_col, 'LineWidth', 1.2, 'LineStyle', '--');
    grid on;  ylabel('u_1(t)');  title(['Control Inputs - ' ttl]);  legend(leg, 'Location', 'best');
    subplot(2,1,2)
    plot(D(1).t, D(1).u2, 'Color', c1_col, 'LineWidth', 1.2, 'LineStyle', '-'); hold on;
    plot(D(2).t, D(2).u2, 'Color', c2_col, 'LineWidth', 1.2, 'LineStyle', '--');
    grid on;  ylabel('u_2(t)');  xlabel('Time [s]');  legend(leg, 'Location', 'best');
    save_png(f, fullfile(opts.outdir, 'control_inputs.png'));

    % ---- norm of z --------------------------------------------------------------------------
    f = figure('Name', 'Norm of z', 'NumberTitle', 'off', 'Visible', vis);
    plot(D(1).t, D(1).nz, 'Color', c1_col, 'LineWidth', 1.5, 'LineStyle', '-'); hold on;
    plot(D(2).t, D(2).nz, 'Color', c2_col, 'LineWidth', 1.5, 'LineStyle', '--');
    grid on;  xlabel('Time [s]');  ylabel('||z(t)|| / ||z_0||');
    title(['Norm of derived z(t) - ' ttl]);  legend(leg, 'Location', 'best');
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
