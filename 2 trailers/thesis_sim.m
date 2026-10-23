function R = thesis_sim(task, opts)
%THESIS_SIM  Data and figures of thesis Chapter 7 (simulation), TWO trailers
%   (r = 4). Twin of "1 trailer/thesis_sim.m".
%
%   thesis_sim                  all tasks, with figures (written to P.fig_dir)
%   thesis_sim('design')        the design of Chapter 6 on the approximation
%                               and on the exact kinematics, side by side
%   thesis_sim('validity')      how far the approximation can be trusted as
%                               the period epsilon grows
%   thesis_sim('adaptive')      fixed gain against the gain recomputed every
%                               interval (approximation)
%   thesis_sim('envelope')      which start poses the exact system converges from
%   thesis_sim('stepstab')      contraction factor per interval against eps*gamma (limits ignored)
%   thesis_sim('long')          (slow; opts.epsilon = other period, not in 'all') interval by interval on the exact
%                               kinematics up to opts.t_end, against the gradient-flow estimate
%   R = thesis_sim(task, opts)  opts.plots = false  : no figures, no files
%                               opts.eps_list       : periods for 'validity'
%                               opts.scales         : start scales for 'envelope'
%                               opts.t_end          : run length for 'design'
%
%   All simulations use the functions of trailer2.m (simulate_case,
%   simulate_case_physical, simulate_case_adaptive) through
%   trailer2('name', ...). The design (period, multipliers, Lyapunov
%   function) is the same block as in thesis_figures.m; the gain is the
%   largest one admissible along the ray (Chapter 6, eq. (6.3)).
%
%   The output grid is dt = min(0.02 s, 1/40 of the fastest oscillation);
%   u1 attains its peak at the start of every interval, which is on the
%   grid, so the peak of |u1| is exact; the peak of |u2| is accurate to a
%   few percent.

    if nargin < 1 || isempty(task), task = 'all'; end
    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'plots'), opts.plots = true; end
    % phase offsets of the three oscillating groups (test of the phase choice)
    if isfield(opts, 'phase'), setappdata(0, 'ctrl_phase', opts.phase(:).');
    else, setappdata(0, 'ctrl_phase', [0 0 0]); end
    cleanup_phase = onCleanup(@() setappdata(0, 'ctrl_phase', [0 0 0]));
    P  = params();
    Tr = @(name, varargin) trailer2(name, varargin{:});

    % ---------------- the design (same block as thesis_figures.m) ----------
    DES.epsilon = 3;  DES.k = [14 5 1];  DES.delta = 0.005;
    h_use     = @(z) 1/sqrt(norm(z)^2 + DES.delta^2);
    lyap_name = sprintf('sqrt, delta = %g', DES.delta);

    % ---------------- model adapter (the ONLY part that differs from the
    % two-trailer file) ---------------------------------------------------
    MODEL = '2tr';  LABEL = 'two trailers (r = 4)';  NANG = 3;
    geo = struct('d0', P.d0, 'l1', P.l1, 'd1', P.d1, 'l2', P.l2);
    [alpha, kappa, theta] = offhooked_coeffs(geo.d0, geo.d1, geo.l1, geo.l2);
    assert(norm([alpha kappa theta] - [P.alpha P.kappa P.theta]) < 1e-14, 'coefficients differ from params.m');
    x0 = P.x0(:);                                                  % radians here
    x0at = @(s, sg) [sg*s*0.2; sg*s*0.2; s*deg2rad(45); s*deg2rad(45); s*deg2rad(45)];
    z0of = @(x) x_to_z(x, geo.d0, geo.d1, geo.l1, geo.l2);
    simA = @(z0, gam, ep, T, dt) Tr('simulate_case', z0, alpha, kappa, theta, gam, ep, 2*pi/ep, ...
                                    DES.k(1), DES.k(2), DES.k(3), 0, T, dt, h_use);
    simE = @(x, gam, ep, T, dt) Tr('simulate_case_physical', x, geo.d0, geo.d1, geo.l1, geo.l2, ...
                                   alpha, kappa, theta, gam, ep, 2*pi/ep, ...
                                   DES.k(1), DES.k(2), DES.k(3), 0, T, dt, h_use);
    simD = @(z0, ep, T, dt, U1, U2, cap) Tr('simulate_case_adaptive', z0, alpha, kappa, theta, ep, ...
                                   2*pi/ep, DES.k(1), DES.k(2), DES.k(3), 0, T, dt, h_use, [], [], U1, U2, cap);
    gmax = @(z, om, u1, u2, h) Tr('gamma_max', z, alpha, kappa, theta, om, ...
                                  DES.k(1), DES.k(2), DES.k(3), u1, u2, h);

    U1 = P.u1_max;  U2 = P.u2_max;  cap = P.cap_product;
    z0 = z0of(x0);  r0 = norm(z0);
    assert(norm(z0 - P.z0) < 1e-12, 'z0 differs from params.m');
    gam_of = @(ep) ray_gain(gmax, z0, ep, U1, U2, cap, h_use);
    dt_of  = @(ep) min(0.02, 2*pi/(max(DES.k)*2*pi/ep)/40);

    R = struct();
    R.model = MODEL;  R.label = LABEL;  R.des = DES;  R.z0 = z0;  R.x0 = x0;
    R.u1max = U1;  R.u2max = U2;
    outdir = P.fig_dir;
    if opts.plots && ~exist(outdir, 'dir'), mkdir(outdir); end
    run = @(name) any(strcmp(task, {'all', name}));

    fprintf('\n=== thesis_sim: %s, design eps = %g s, k = (%s), V = %s ===\n', ...
            LABEL, DES.epsilon, num2str(DES.k), lyap_name);

    if run('design')
        ep = DES.epsilon;  gam = gam_of(ep);  dt = dt_of(ep);
        if isfield(opts, 't_end'), T = opts.t_end; else, T = 60; end
        fprintf('\n-- design run: gamma = %.4e, dt = %.4f s, t_end = %g s\n', gam, dt, T);
        [ta, za, u1a, u2a]       = simA(z0, gam, ep, T, dt);
        [te, xe, ze, u1e, u2e]   = simE(x0, gam, ep, T, dt);
        D = struct('gamma', gam, 'dt', dt, 'T', T, 'eps', ep, ...
                   'ta', ta, 'za', za, 'u1a', u1a, 'u2a', u2a, ...
                   'te', te, 'xe', xe, 'ze', ze, 'u1e', u1e, 'u2e', u2e);
        na = vecnorm(za, 2, 2);  ne = vecnorm(ze, 2, 2);
        lv = [0.1 0.01 0.001];
        D.t_cross_a = arrayfun(@(q) t_cross(ta, na, q*r0), lv);
        D.t_cross_e = arrayfun(@(q) t_cross(te, ne, q*r0), lv);
        n = min(numel(ta), numel(te));
        dev = vecnorm(za(1:n,:) - ze(1:n,:), 2, 2);
        [D.dev_max, i] = max(dev);  D.dev_rel = D.dev_max / r0;  D.dev_at = ta(i);
        D.u1_max_a = max(abs(u1a));  D.u2_max_a = max(abs(u2a));
        D.u1_max_e = max(abs(u1e));  D.u2_max_e = max(abs(u2e));
        D.final_a = norm(za(end,:));  D.final_e = norm(ze(end,:));
        D.theta_final_deg = rad2deg(xe(end, 4:5));
        D.theta_max_deg   = rad2deg(max(abs(xe(:, 4:5))));
        fprintf('  ||z|| to 10%% / 1%% / 0.1%% of ||z0||  approx %s s   exact %s s\n', ...
                mat2str(D.t_cross_a, 4), mat2str(D.t_cross_e, 4));
        fprintf('  largest ||z_exact - z_approx|| = %.4f (%.1f %% of ||z0||) at t = %.1f s\n', ...
                D.dev_max, 100*D.dev_rel, D.dev_at);
        fprintf('  max |u1| approx %.4f exact %.4f (limit %.2f);  max |u2| approx %.4f exact %.4f (limit %.2f)\n', ...
                D.u1_max_a, D.u1_max_e, U1, D.u2_max_a, D.u2_max_e, U2);
        fprintf('  final ||z|| approx %.3e exact %.3e;  final theta1, theta2 %s deg, largest |theta| %s deg\n', ...
                D.final_a, D.final_e, mat2str(D.theta_final_deg, 4), mat2str(D.theta_max_deg, 4));
        R.design = D;
        if opts.plots, plot_design(D, MODEL, outdir, U1, U2, r0); end
    end

    if run('validity')
        if isfield(opts, 'eps_list'), EL = opts.eps_list; else, EL = [2 3 4 6 8 26]; end
        N = 12;
        V = struct('eps', EL, 'gamma', zeros(size(EL)), 'e_rel', zeros(size(EL)), ...
                   'prog_a', zeros(size(EL)), 'prog_e', zeros(size(EL)), ...
                   'u1_e', zeros(size(EL)), 'u2_e', zeros(size(EL)), ...
                   'u1_a', zeros(size(EL)), 'ang_max', zeros(size(EL)), ...
                   'excursion', zeros(size(EL)));
        fprintf('\n-- validity: %d intervals per period, gain = largest admissible\n', N);
        fprintf('  eps    gamma       dev/||z0||  ||z(T)||/||z0||  approx  exact   |u1|e/lim  |u2|e/lim  max|angle| [deg]  excursion [m]\n');
        for i = 1:numel(EL)
            ep = EL(i);  gam = gam_of(ep);  dt = dt_of(ep);  T = N*ep;
            [ta, za, u1a]          = simA(z0, gam, ep, T, dt);
            [te, xe, ze, u1e, u2e] = simE(x0, gam, ep, T, dt);
            n = min(numel(ta), numel(te));
            V.gamma(i)  = gam;
            V.e_rel(i)  = max(vecnorm(za(1:n,:) - ze(1:n,:), 2, 2)) / r0;
            V.prog_a(i) = norm(za(end,:)) / r0;
            V.prog_e(i) = norm(ze(end,:)) / r0;
            V.u1_e(i)   = max(abs(u1e)) / U1;  V.u2_e(i) = max(abs(u2e)) / U2;
            V.u1_a(i)   = max(abs(u1a)) / U1;
            V.ang_max(i) = rad2deg(max(max(abs(xe(:, 3:end)))));
            % largest displacement of the robot within ONE interval, over all intervals
            ex = 0;
            for j = 1:N
                sel = te >= (j-1)*ep - 1e-9 & te <= j*ep + 1e-9;
                p = xe(sel, 1:2);
                ex = max(ex, max(vecnorm(p - p(1,:), 2, 2)));
            end
            V.excursion(i) = ex;
            fprintf('  %5.1f  %.3e   %8.3f     %8.3f   %8.3f     %6.3f    %6.3f     %8.1f        %.3f\n', ...
                    ep, gam, V.e_rel(i), V.prog_a(i), V.prog_e(i), V.u1_e(i), V.u2_e(i), V.ang_max(i), ex);
        end
        R.validity = V;
        if opts.plots, plot_validity(V, MODEL, outdir); end
    end

    if run('adaptive')
        ep = DES.epsilon;  gam = gam_of(ep);  dt = dt_of(ep);
        if isfield(opts, 't_end'), T = opts.t_end; else, T = 60; end
        fprintf('\n-- fixed against adaptive gain (approximation), t_end = %g s\n', T);
        [tf, zf, u1f, u2f]            = simA(z0, gam, ep, T, dt);
        [tq, zq, u1q, u2q, gq]        = simD(z0, ep, T, dt, U1, U2, cap);
        A = struct('gamma_fixed', gam, 'tf', tf, 'zf', zf, 'tq', tq, 'zq', zq, 'gq', gq, ...
                   'u1f', u1f, 'u2f', u2f, 'u1q', u1q, 'u2q', u2q);
        lv = [0.1 0.01 0.001];
        A.t_cross_f = arrayfun(@(q) t_cross(tf, vecnorm(zf,2,2), q*r0), lv);
        A.t_cross_q = arrayfun(@(q) t_cross(tq, vecnorm(zq,2,2), q*r0), lv);
        fprintf('  fixed    : to 10/1/0.1 %%: %s s, max |u1| %.4f, max |u2| %.4f\n', mat2str(A.t_cross_f,4), max(abs(u1f)), max(abs(u2f)));
        fprintf('  adaptive : to 10/1/0.1 %%: %s s, max |u1| %.4f, max |u2| %.4f\n', mat2str(A.t_cross_q,4), max(abs(u1q)), max(abs(u2q)));
        fprintf('  adaptive gamma: min %.3e  max %.3e  (fixed %.3e);  final ||z|| fixed %.2e adaptive %.2e\n', ...
                min(gq), max(gq), gam, norm(zf(end,:)), norm(zq(end,:)));
        R.adaptive = A;
        if opts.plots, plot_adaptive(A, MODEL, outdir, r0); end
    end

    if run('envelope')
        if isfield(opts, 'scales'), SC = opts.scales; else, SC = [1 0.75 0.5 0.4 0.3 0.2]; end
        ep = DES.epsilon;  gam = gam_of(ep);  dt = dt_of(ep);
        fprintf('\n-- envelope of start poses x0 = s*(+-0.2 m, +-0.2 m, 45 deg, 45 deg), exact kinematics, 60 s\n');
        fprintf('  sign  s     |x0|[m]   ||z(60 s)||/||z0||   max|u1|/lim   theta1 end [deg]\n');
        E = struct('sign', [], 's', [], 't001', [], 'u1', [], 'theta', []);
        for sg = [1 -1]
            for s = SC
                if sg < 0 && s ~= 1, continue, end
                x = x0at(s, sg);  zz = z0of(x);
                [te, xe, ze, u1e] = simE(x, gam, ep, 60, dt);
                tt = norm(ze(end,:)) / norm(zz);
                E.sign(end+1) = sg;  E.s(end+1) = s;  E.t001(end+1) = tt;
                E.u1(end+1) = max(abs(u1e))/U1;  E.theta(end+1) = rad2deg(xe(end,4));
                fprintf('  %+d    %.2f   %.3f     %8.1f       %6.3f        %8.1f\n', sg, s, ...
                        norm(x(1:2)), tt, E.u1(end), E.theta(end));
            end
        end
        R.envelope = E;
    end

    % step stability of the sampled scheme (approximation, quadratic V, h = 1):
    % the contraction factor per interval is 1 - eps*gamma*h, so the scheme
    % converges for 0 < eps*gamma*h < 2 and is fastest at 1. Amplitude limits
    % are NOT applied here (this is the contraction condition alone).
    if run('stepstab')
        ep = DES.epsilon;  N = 10;  PL = [0.4 0.8 1.0 1.2 1.4 1.6 1.9 2.1 2.6];
        S = struct('prod', PL, 'ratio', zeros(size(PL)), 'theory', abs(1 - PL).^N);
        fprintf('\n-- step stability: eps = %g s, %d intervals, quadratic V (h = 1), limits ignored\n', ep, N);
        fprintf('  eps*gamma   ||z_N||/||z_0||   (1-eps*gamma)^N\n');
        for i = 1:numel(PL)
            gam = PL(i)/ep;
            [~, zq] = Tr('simulate_case', z0, alpha, kappa, theta, gam, ep, 2*pi/ep, DES.k(1), DES.k(2), DES.k(3), 0, N*ep, 0.02, @(z) 1);
            S.ratio(i) = norm(zq(end,:)) / r0;
            fprintf('  %7.2f    %12.3e     %12.3e\n', PL(i), S.ratio(i), S.theory(i));
        end
        R.stepstab = S;
    end

    % 'long' is NOT part of 'all' (it takes long): the exact kinematics, one
    % interval at a time, only the state at the end of each interval is kept
    % (the amplitudes are frozen at the start of an interval and the periods
    % are whole multiples of the oscillations, so restarting every interval
    % changes nothing). Compared with the gradient-flow estimate of Chapter 6.
    if strcmp(task, 'long')
        ep = DES.epsilon;  if isfield(opts, 'epsilon'), ep = opts.epsilon; end
        gam = gam_of(ep);  dt = dt_of(ep);
        if isfield(opts, 'scale') && opts.scale ~= 1
            sc = opts.scale;                      % smaller start: x0 -> sc*x0 (positions and angles)
            x0 = sc*x0;  z0 = z0of(x0);  r0 = norm(z0);
            gam = ray_gain(gmax, z0, ep, U1, U2, cap, h_use);   % gain tuned along the ray of THIS start
            fprintf('\n-- start scaled by %g: ||z0|| = %.4f\n', sc, r0);
        end
        if isfield(opts, 'gain_factor')               % run with a fraction of the admissible gain (test)
            gam = gam*opts.gain_factor;
            fprintf('\n-- gain multiplied by %g\n', opts.gain_factor);
        end
        if isfield(opts, 't_end'), T = opts.t_end; else, T = 12*ep; end
        N = ceil(T/ep - 1e-9);
        fprintf('\n-- long run: %d intervals of %g s (t_end = %g s), gamma = %.4e\n', N, ep, N*ep, gam);
        xc = x0;  zc = z0;
        tt = (0:N)*ep;  ne = zeros(1, N+1);  na = ne;  u1m = zeros(1, N);
        u1a_m = zeros(1, N);
        angm = 0;
        hitm = [0 0];  hit_over = 0;  nsamp = 0;     % largest hitch angles over ALL samples, samples above 1.40 rad
        ne(1) = 1;  na(1) = 1;
        t0 = tic;
        for j = 1:N
            [~, xs, ~, u1s] = simE(xc, gam, ep, ep, dt);
            xc = xs(end, :).';  ne(j+1) = norm(z0of(xc)) / r0;  u1m(j) = max(abs(u1s));
            angm = max(angm, rad2deg(max(abs(xc(3:end)))));
            hh = abs([xs(:,3) - xs(:,4), xs(:,4) - xs(:,5)]);    % car-trailer1, trailer1-trailer2 [rad]
            hitm = max(hitm, max(hh, [], 1));
            hit_over = hit_over + sum(any(hh > 1.40, 2));  nsamp = nsamp + size(hh, 1);
            if ~isfield(opts, 'approx') || opts.approx
                [~, zs, u1as] = simA(zc, gam, ep, ep, dt);
                u1a_m(j) = max(abs(u1as));
                zc = zs(end, :).';  na(j+1) = norm(zc) / r0;
            end
            if mod(j, max(1, round(N/10))) == 0
                fprintf('   %5.1f %%  t = %9.0f s  exact %.4f  approx %.4f   (%.0f s elapsed)\n', ...
                        100*j/N, tt(j+1), ne(j+1), na(j+1), toc(t0));
            end
        end
        fr = logspace(0, log10(max(min([ne na]), 1e-3)), 80);
        te = Tr('eta_ray', z0, gam, h_use, fr);
        est = exp(interp1(te, log(fr), tt, 'linear', 'extrap'));
        L = struct('t', tt, 'exact', ne, 'approx', na, 'est', est, 'u1max', max(u1m)/U1, ...
                   'gamma', gam, 'eps', ep, 'ang_end_max', angm, 'u1_int', u1m/U1, 'u1_int_approx', u1a_m/U1);
        L.err_est = max(abs(ne - est));      % absolute, in units of ||z0||
        L.err_ap  = max(abs(ne - na));
        fprintf('  at t = %g s: exact %.5f, approx %.5f, gradient-flow estimate %.5f\n', ...
                tt(end), ne(end), na(end), est(end));
        fprintf('  largest difference (in %% of ||z0||) exact/estimate %.2f, exact/approx %.2f; max |u1| = %.4f of the limit\n', ...
                100*L.err_est, 100*L.err_ap, L.u1max);
        [mx, jx] = max(u1m/U1);
        fprintf('  peak |u1| of the exact run is %.4f of the limit, in interval %d (t = %g s);  approximation: %.4f\n', ...
                mx, jx, (jx-1)*ep, max(u1a_m)/U1);
        fprintf('  intervals with |u1| above the limit (exact): %d of %d\n', sum(u1m/U1 > 1 + 1e-9), N);
        fprintf('  largest |angle| at the end of an interval: %.1f deg\n', angm);
        fprintf('  largest hitch angle over all samples: car-trailer1 %.3f rad, trailer1-trailer2 %.3f rad; samples above 1.40 rad: %.1f %%\n', ...
                hitm(1), hitm(2), 100*hit_over/max(nsamp,1));
        fprintf('  exact model reaches 5 %% / 1 %% / 0.1 %% of ||z0|| at t = %s s\n', ...
                mat2str(arrayfun(@(q) t_cross(tt, ne, q), [0.05 0.01 0.001]), 5));
        R.long = L;
        if opts.plots, plot_long(L, MODEL, outdir); end
    end

    if opts.plots
        matname = ['sim_' MODEL '.mat'];
        if isfield(R, 'long') && ~isfield(R, 'design')      % a lone 'long' run must not overwrite the full .mat
            matname = ['sim_' MODEL '_long.mat'];
        end
        save(fullfile(outdir, matname), '-struct', 'R');
        fprintf('\nfigures and sim_%s.mat written to %s\n', MODEL, outdir);
    end
end

% ---- helpers ----

function g = ray_gain(gmax, z0, epsilon, U1, U2, cap, h_fun)
%RAY_GAIN  Largest gain admissible along the whole ray z0 -> 0
%   (amplitude at every point, and eps*gamma*h <= cap); as in thesis_figures.m
    r0 = norm(z0);  zh = z0(:)/r0;  om = 2*pi/epsilon;
    g = inf;  hmax = 0;
    for r = r0*logspace(-3, 0, 200)
        g    = min(g, gmax(r*zh, om, U1, U2, h_fun));
        hmax = max(hmax, h_fun(r*zh));
    end
    g = min(g, cap/(epsilon*hmax));
end

function t = t_cross(t, nz, level)
%T_CROSS  first time at which the norm nz drops to level (NaN if never)
    i = find(nz <= level, 1, 'first');
    if isempty(i), t = NaN; return, end
    if i == 1, t = t(1); return, end
    f = (nz(i-1) - level) / (nz(i-1) - nz(i));
    t = t(i-1) + f*(t(i) - t(i-1));
end

function C = palette()
    C.a = [0.00 0.45 0.74];  C.e = [0.85 0.33 0.10];  C.g = [0.45 0.45 0.45];
end

function f = new_fig()
    f = figure('Color', 'w');
    set(f, 'DefaultAxesFontSize', 9, 'DefaultTextFontSize', 9);
end

function save_fig(f, name, outdir, w_cm, h_cm)
    set(f, 'Units', 'centimeters', 'Position', [2 2 w_cm h_cm], ...
           'PaperUnits', 'centimeters', 'PaperSize', [w_cm h_cm], ...
           'PaperPosition', [0 0 w_cm h_cm]);
    try, set(findall(f, 'Type', 'axes'), 'Toolbar', []); catch, end   % no axes toolbar in the PDF
    print(f, fullfile(outdir, name), '-dpdf', '-painters');
    close(f);
end

function plot_design(D, MODEL, outdir, U1, U2, r0)
    C = palette();
    if D.T >= 7200, u = 3600; ul = 'h'; elseif D.T >= 600, u = 60; ul = 'min'; else, u = 1; ul = 's'; end
    % convergence
    f = new_fig();
    semilogy(D.ta/u, max(vecnorm(D.za,2,2)/r0, 1e-12), 'Color', C.a, 'LineWidth', 1.6); hold on
    semilogy(D.te/u, max(vecnorm(D.ze,2,2)/r0, 1e-12), '--', 'Color', C.e, 'LineWidth', 1.6);
    grid on; xlabel(['t  [' ul ']']); ylabel('||z(t)|| / ||z_0||');
    legend('approximation', 'exact kinematics', 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
    save_fig(f, ['sim_' MODEL '_conv'], outdir, 12, 7);
    % controls over the first three intervals
    f = new_fig();
    m = D.ta <= 3*D.eps + 1e-9;  me = D.te <= 3*D.eps + 1e-9;
    subplot(2,1,1)
    plot(D.ta(m), D.u1a(m), 'Color', C.a, 'LineWidth', 1.2); hold on
    plot(D.te(me), D.u1e(me), '--', 'Color', C.e, 'LineWidth', 1.2);
    hline(U1, C.g); hline(-U1, C.g);
    grid on; ylabel('u_1  [m/s]');
    legend('approximation', 'exact kinematics', 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
    subplot(2,1,2)
    plot(D.ta(m), D.u2a(m), 'Color', C.a, 'LineWidth', 1.2); hold on
    plot(D.te(me), D.u2e(me), '--', 'Color', C.e, 'LineWidth', 1.2);
    hline(U2, C.g); hline(-U2, C.g);
    grid on; xlabel('t  [s]'); ylabel('u_2  [rad/s]');
    save_fig(f, ['sim_' MODEL '_ctrl'], outdir, 12, 9);
    % path of the robot (exact)
    f = new_fig();
    plot(D.xe(:,1), D.xe(:,2), 'Color', C.e, 'LineWidth', 1.2); hold on
    plot(D.xe(1,1), D.xe(1,2), 'o', 'Color', C.e, 'MarkerFaceColor', C.e);
    plot(0, 0, 'k+', 'MarkerSize', 8, 'LineWidth', 1.2);
    grid on; axis equal; xlabel('x_1  [m]'); ylabel('x_2  [m]');
    save_fig(f, ['sim_' MODEL '_path'], outdir, 8, 7);
end

function plot_validity(V, MODEL, outdir)
    C = palette();
    f = new_fig();
    subplot(1,2,1)
    loglog(V.eps, max(V.e_rel, 1e-4), 'o-', 'Color', C.a, 'LineWidth', 1.5);
    grid on; xlabel('\epsilon  [s]'); ylabel('max ||z_{exact} - z_{approx}|| / ||z_0||');
    subplot(1,2,2)
    semilogx(V.eps, V.prog_a, 'o-', 'Color', C.a, 'LineWidth', 1.5); hold on
    semilogx(V.eps, V.prog_e, 's--', 'Color', C.e, 'LineWidth', 1.5);
    grid on; xlabel('\epsilon  [s]'); ylabel('||z(12 \epsilon)|| / ||z_0||');
    legend('approximation', 'exact kinematics', 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
    save_fig(f, ['sim_' MODEL '_validity'], outdir, 14, 6);
end

function plot_adaptive(A, MODEL, outdir, r0)
    C = palette();
    f = new_fig();
    subplot(1,2,1)
    semilogy(A.tf, vecnorm(A.zf,2,2)/r0, 'Color', C.a, 'LineWidth', 1.5); hold on
    semilogy(A.tq, vecnorm(A.zq,2,2)/r0, '--', 'Color', C.e, 'LineWidth', 1.5);
    grid on; xlabel('t  [s]'); ylabel('||z(t)|| / ||z_0||');
    legend('fixed gain', 'adaptive gain', 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
    subplot(1,2,2)
    semilogy(A.tq, A.gq, 'Color', C.e, 'LineWidth', 1.5); hold on
    hline(A.gamma_fixed, C.a);
    grid on; xlabel('t  [s]'); ylabel('\gamma');
    save_fig(f, ['sim_' MODEL '_adaptive'], outdir, 14, 6);
end

function plot_long(L, MODEL, outdir)
    C = palette();
    if L.t(end) >= 7200, u = 3600; ul = 'h'; elseif L.t(end) >= 600, u = 60; ul = 'min'; else, u = 1; ul = 's'; end
    f = new_fig();
    semilogy(L.t/u, L.est, '-', 'Color', C.g, 'LineWidth', 1.4); hold on
    semilogy(L.t/u, L.approx, '-', 'Color', C.a, 'LineWidth', 1.4);
    semilogy(L.t/u, L.exact, '--', 'Color', C.e, 'LineWidth', 1.4);
    grid on; xlabel(['t  [' ul ']']); ylabel('||z|| / ||z_0||  at the end of each interval');
    legend('gradient-flow estimate', 'approximation', 'exact kinematics', 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
    save_fig(f, ['sim_' MODEL '_long'], outdir, 12, 7);
end

function hline(y, col)
    h = yline(y, ':');
    set(h, 'Color', col);
end
