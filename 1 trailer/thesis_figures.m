function D = thesis_figures(do_plots)
%THESIS_FIGURES  Data and figures of thesis Chapter 6 (feasibility limits),
%   ONE trailer (r = 3, four privileged coordinates).
%
%   thesis_figures            compute everything, draw the figures, write the
%                             PDFs into P.fig_dir and the numbers into
%                             feas_1tr.mat next to them
%   D = thesis_figures(false) compute only (no figures, no files)
%
%   Everything the chapter quotes comes out of this one call, so the text
%   and the pictures cannot drift apart. The numbers use the SAME functions
%   as the rest of the study (trailer1.m: gamma_max, control_bound,
%   eta_ray, x_to_z), reached through trailer1('name', ...). Nothing is
%   re-implemented here except the bookkeeping.
%
%   The two-trailer twin of this file is "2 trailers/thesis_figures.m"; the
%   figure that sets r = 2, 3, 4 side by side is made by
%   thesis_figures_compare.m after both have run.
%
%   Conventions
%     * "quadratic V" (h = 1) is used for the sweeps over epsilon, geometry
%       and hardware limits, because there the answer does not depend on
%       the extra parameter delta; the Lyapunov functions themselves are
%       compared in fig "V".
%     * gain = the larger of nothing: effective gain is
%           gamma_eff = min( gamma_amp , cap_product/(epsilon*hmax) )
%       with gamma_amp = largest gain for which control_bound stays inside
%       both velocity limits, and the second term the contraction cap.
%     * time to 10 % of ||z0|| comes from the gradient-flow estimate
%       (eta_ray), valid for any radial V.

    if nargin < 1, do_plots = true; end
    P  = params();
    Tr = @(name, varargin) trailer1(name, varargin{:});

    % ---------------- the design of the experiment ---------------------------
    % Period, multipliers and the Lyapunov function of Chapter 6 are fixed HERE
    % (not read from params.m), so the chapter does not change when params.m is
    % edited for other purposes. Geometry, start pose x0, limits and the cap
    % still come from params.m. The gain is not an input: it is the largest one
    % admissible along the ray (ray_gain), which is how params.m was tuned.
    DES.epsilon = 5;  DES.k = [3 2];     DES.c = 10;
    DES.delta = 0.02;                    % only for the sqrt candidate in the V comparison
    h_use     = @(z) 1/(1 + DES.c*norm(z)^2);
    lyap_name = sprintf('log, c = %g (Victoria)', DES.c);

    % ---------------- model adapter (the ONLY part that differs from the
    % two-trailer file) --------------------------------------------------
    MODEL = '1tr';  R = 3;  LABEL = 'one trailer (r = 3)';
    geo0 = struct('d0', P.d0, 'l1', P.l1, 'd1', 0, 'l2', 0);
    x0   = [P.x0(1); P.x0(2); deg2rad(P.x0(3)); deg2rad(P.x0(4))];  % P.x0 is in degrees
    kvec = DES.k;
    coef = @(g) g.d0*(g.d0 + g.l1)/2;                       % alpha (params.m, one_trailer_coeff)
    z0of = @(g) Tr('x_to_z', x0, g.d0, g.l1);
    z0sc = @(g, s) Tr('x_to_z', [s*x0(1:2); x0(3:end)], g.d0, g.l1);   % start distance scaled by s
    gmaxk = @(z, g, omega, k, u1, u2, h) Tr('gamma_max', z, coef(g), omega, ...
                                            k(1), k(2), u1, u2, h);
    cbnd = @(z, g, gam, omega, h) cb_wrap(Tr, z, coef(g), gam, omega, kvec, h);
    KTAB = [3 2; 5 2; 7 2; 2 3];                       % multipliers compared in section 7 (first row = in use)
    gmax = @(z, g, omega, u1, u2, h) gmaxk(z, g, omega, kvec, u1, u2, h);

    TU = 60;    TUL = 'min';
    U1 = P.u1_max;  U2 = P.u2_max;  cap = P.cap_product;
    eps0 = DES.epsilon;  hq = @(z) 1;
    z0 = P.z0;  r0 = norm(z0);

    % the adapter must reproduce the geometry-dependent pieces of params.m
    assert(abs(coef(geo0) - P.alpha) < 1e-14, 'alpha in adapter differs from params.m');
    assert(norm(z0of(geo0) - z0) < 1e-12,     'z0 in adapter differs from params.m');

    D.model = MODEL;  D.r = R;  D.label = LABEL;
    D.geo = geo0;  D.z0 = z0;  D.k = kvec;  D.eps_use = eps0;
    D.lyap_use = lyap_name;  D.delta = DES.delta;
    D.u1max = U1;  D.u2max = U2;  D.cap = cap;

    % ---------------- 1. epsilon sweep, quadratic V ---------------------------
    eps_list = [0.25 0.5 1 2 3 4 5 6 8 12 16 24 32];
    n = numel(eps_list);
    g_amp = zeros(1, n);
    for i = 1:n
        om = 2*pi/eps_list(i);
        g_amp(i) = gmax(z0, geo0, om, U1, U2, hq);
    end
    g_cap = cap ./ eps_list;
    g_eff = min(g_amp, g_cap);
    D.eps_list = eps_list;  D.g_amp = g_amp;  D.g_cap = g_cap;  D.g_eff = g_eff;
    D.t10 = log(10) ./ g_eff;                % [s], gradient flow, h = 1
    D.t10_unit = 's';

    % is z0 really the worst point of the ray for the quadratic V?
    % (it is what the shortcut above assumes)
    rs = r0 * logspace(-3, 0, 60);  worst = inf;
    om0 = 2*pi/eps0;
    for r = rs
        worst = min(worst, gmax(r*z0/r0, geo0, om0, U1, U2, hq));
    end
    D.ray_check = worst / gmax(z0, geo0, om0, U1, U2, hq);   % must be 1

    % ---------------- 2. which terms use the budget (in-use design) -----------
    D.gamma_use = ray_gain(gmax, z0, geo0, eps0, U1, U2, cap, h_use);
    [b1, b2, t1, t2] = cbnd(z0, geo0, D.gamma_use, 2*pi/eps0, h_use);
    D.terms_u1 = t1 / U1;   D.terms_u2 = t2 / U2;
    D.terms_names = {'a (constant)', 'T12', 'T112'};
    D.bound_u1 = b1 / U1;   D.bound_u2 = b2 / U2;

    % ---------------- 3. Lyapunov functions ----------------------------------
    CAND = { 'quadratic',                @(z) 1; ...
             'log, c = 10',              @(z) 1/(1 + 10*norm(z)^2); ...
             'log, c = 100',             @(z) 1/(1 + 100*norm(z)^2); ...
             ['sqrt, \delta = ' num2str(DES.delta)], @(z) 1/sqrt(norm(z)^2 + DES.delta^2) };
    nC = size(CAND, 1);
    fr = logspace(0, -3, 80);
    D.V_name = CAND(:,1);  D.V_gamma = zeros(nC,1);  D.V_t = zeros(nC, numel(fr));
    D.V_frac = fr;
    for c = 1:nC
        h = CAND{c,2};
        g = ray_gain(gmax, z0, geo0, eps0, U1, U2, cap, h);
        D.V_gamma(c) = g;
        D.V_t(c,:) = Tr('eta_ray', z0, g, h, fr);
    end
    D.V_t10  = interp_t(D.V_t, fr, 0.1);
    D.V_t001 = interp_t(D.V_t, fr, 0.001);

    % ---------------- 4. geometry sweeps (quadratic V, epsilon in use) ---------
    d0_list = [0.02 0.04 0.06 0.08 0.10 0.115 0.14 0.17 0.20 0.25 0.30];
    ls_list = [0.5 0.75 1 1.5 2 3 4];
    D.d0_list = d0_list;  D.ls_list = ls_list;
    [D.d0_g, D.d0_t10, D.d0_ampfrac] = geo_sweep(gmax, z0of, geo0, 'd0', d0_list, ...
                                       eps0, U1, U2, cap);
    [D.ls_g, D.ls_t10, D.ls_ampfrac] = geo_sweep(gmax, z0of, geo0, 'l',  ls_list, ...
                                       eps0, U1, U2, cap);

    % ---------------- 5. hardware limits (quadratic V, epsilon in use) ---------
    s_list = [0.25 0.5 0.75 1 1.5 2 3 4 6 8];
    D.s_list = s_list;
    D.hw_u1 = zeros(size(s_list));  D.hw_u2 = D.hw_u1;
    for i = 1:numel(s_list)
        ga = gmax(z0, geo0, 2*pi/eps0, s_list(i)*U1, U2, hq);
        D.hw_u1(i) = log(10) / min(ga, cap/eps0);
        ga = gmax(z0, geo0, 2*pi/eps0, U1, s_list(i)*U2, hq);
        D.hw_u2(i) = log(10) / min(ga, cap/eps0);
    end

    % ---------------- 7. multipliers k (quadratic V, epsilon in use) ------------
    D.ktab = KTAB;  D.k_g = zeros(size(KTAB,1), 1);  D.k_gray = D.k_g;  D.k_t001 = D.k_g;
    for i = 1:size(KTAB,1)
        D.k_g(i) = gmaxk(z0, geo0, 2*pi/eps0, KTAB(i,:), U1, U2, hq);
        % same multipliers with the V in use, along the whole ray (what params.m is tuned by)
        gk = @(z, g, om, u1, u2, h) gmaxk(z, g, om, KTAB(i,:), u1, u2, h);
        D.k_gray(i,1) = ray_gain(gk, z0, geo0, eps0, U1, U2, cap, h_use);
        D.k_t001(i,1) = Tr('eta_ray', z0, D.k_gray(i), h_use, 0.001);
    end

    % ---------------- 8. start distance (quadratic V, epsilon in use) ----------
    % position part of x0 scaled by s, the angles are kept
    D.sc_list = [0.25 0.5 1 2 4 8];
    D.sc_g = zeros(size(D.sc_list));  D.sc_z0 = D.sc_g;
    for i = 1:numel(D.sc_list)
        zs = z0sc(geo0, D.sc_list(i));
        D.sc_z0(i) = norm(zs);
        D.sc_g(i)  = min(gmax(zs, geo0, 2*pi/eps0, U1, U2, hq), cap/eps0);
    end
    D.sc_t10 = log(10) ./ D.sc_g;

    % ---------------- 9. best possible radial design ---------------------------
    % gain chosen at every radius as large as the limits and the cap allow
    % (gamma_adm(r) = min(gamma_amp(r), cap/eps), quadratic h = 1); then
    % dr/dt = -gamma_adm(r) r and t(q) = int dr / (gamma_adm r).
    rr = r0*logspace(-4, 0, 400);  gadm = zeros(size(rr));
    for i = 1:numel(rr)
        gadm(i) = min(gmax(rr(i)*z0/r0, geo0, 2*pi/eps0, U1, U2, hq), cap/eps0);
    end
    D.ad_r = rr;  D.ad_gamma = gadm;
    D.ad_t = zeros(size(fr));
    for i = 1:numel(fr)
        sel = rr >= fr(i)*r0;
        lr = log(rr(sel));  ig = 1 ./ gadm(sel);          % int dr/(g r) = int dln(r) / g
        if numel(lr) < 2
            D.ad_t(i) = 0;
        else
            D.ad_t(i) = sum(diff(lr) .* (ig(1:end-1) + ig(2:end)) / 2);
        end
    end
    D.ad_t10  = interp1(log(fr), D.ad_t, log(0.1));
    D.ad_t001 = interp1(log(fr), D.ad_t, log(0.001));

    % ---------------- 6. printout ---------------------------------------------
    fprintf('\n=== thesis_figures: %s ===\n', LABEL);
    fprintf('geometry d0=%.3f l1=%.3f, k=(%s), eps in use = %g s, ||z0|| = %.4f\n', ...
            geo0.d0, geo0.l1, num2str(kvec), eps0, r0);
    fprintf('ray check (worst point of the ray / z0, quadratic V; 1 = z0 is the worst): %.6f\n', D.ray_check);
    fprintf('\n  eps [s]   gamma_amp    gamma_cap    t(10%%) [min]\n');
    for i = 1:n
        fprintf('  %6.2f   %10.3e   %10.3e   %10.2f\n', eps_list(i), g_amp(i), g_cap(i), D.t10(i)/60);
    end
    fprintf('\nin-use design: gamma = %.4e (largest admissible along the ray), V = %s\n', D.gamma_use, lyap_name);
    fprintf('terms at z0, in-use design (share of limit):\n');
    fprintf('  u1: %s   sum %.3f\n', mat2str(D.terms_u1, 3), D.bound_u1);
    fprintf('  u2: %s   sum %.3f\n', mat2str(D.terms_u2, 3), D.bound_u2);
    fprintf('\nLyapunov functions at eps = %g:\n', eps0);
    for c = 1:nC
        fprintf('  %-22s gamma %.3e   t(10%%) %8.2f min   t(0.1%%) %8.2f min\n', ...
                strrep(CAND{c,1}, '\delta', 'd'), D.V_gamma(c), D.V_t10(c)/60, D.V_t001(c)/60);
    end
    fprintf('\nd0 sweep (l1 = %.3f fixed): t(10%%) [min], 1 = amplitude-limited\n', geo0.l1);
    for i = 1:numel(d0_list)
        fprintf('  d0 = %.3f   %9.2f   %5.2f\n', d0_list(i), D.d0_t10(i)/60, D.d0_ampfrac(i));
    end
    fprintf('\nl scale sweep (d0 = %.3f fixed): t(10%%) [min]\n', geo0.d0);
    for i = 1:numel(ls_list)
        fprintf('  l x %.2f    %9.2f   %5.2f\n', ls_list(i), D.ls_t10(i)/60, D.ls_ampfrac(i));
    end
    fprintf('\nhardware limits scaled (eps = %g): t(10%%) [min]\n  scale   u1 scaled   u2 scaled\n', eps0);
    for i = 1:numel(s_list)
        fprintf('  %5.2f   %9.2f   %9.2f\n', s_list(i), D.hw_u1(i)/60, D.hw_u2(i)/60);
    end
    fprintf('\nmultipliers k at eps = %g, quadratic V: gamma_amp\n', eps0);
    for i = 1:size(KTAB,1)
        fprintf('  k = (%s)   %.3e   (%.2f x the first row);  V in use along the ray: gamma %.3e, t(0.1%%) %.2f %s\n', ...
                num2str(KTAB(i,:)), D.k_g(i), D.k_g(i)/D.k_g(1), D.k_gray(i), D.k_t001(i)/TU, TUL);
    end
    fprintf('\nstart distance scaled (eps = %g, quadratic V): ||z0||, gamma, t(10%%)\n', eps0);
    for i = 1:numel(D.sc_list)
        fprintf('  x %5.2f   ||z0|| %7.4f   gamma %.3e   t %8.2f %s\n', D.sc_list(i), D.sc_z0(i), D.sc_g(i), D.sc_t10(i)/TU, TUL);
    end
    fprintf('\nbest radial design (gain at every radius as large as allowed, quadratic h): t(10%%) %.2f %s, t(0.1%%) %.2f %s\n', ...
            D.ad_t10/TU, TUL, D.ad_t001/TU, TUL);
    fprintf('  (log V, c = 10 in use: t(10%%) %.2f %s, t(0.1%%) %.2f %s)\n', D.V_t10(2)/TU, TUL, D.V_t001(2)/TU, TUL);
    fprintf('\nquadratic V, t(0.1%%) for every multiplier set at eps = %g:\n', eps0);
    for i = 1:size(KTAB,1)
        fprintf('  k = (%s)   %.2f %s\n', num2str(KTAB(i,:)), log(1000)/D.k_g(i)/TU, TUL);
    end

    if ~do_plots, return, end

    % ---------------- files ---------------------------------------------------
    outdir = P.fig_dir;
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    save(fullfile(outdir, ['feas_' MODEL '.mat']), '-struct', 'D');

    C = palette();
    unit = 60;  ulab = 'min';

    % -- V : convergence for the candidate Lyapunov functions
    f = new_fig();
    for c = 1:size(CAND,1)
        semilogy(D.V_t(c,:)/unit, fr, 'LineWidth', 1.6, 'Color', C.cat(c,:)); hold on
    end
    grid on; xlabel(['t  [' ulab ']']); ylabel('||z(t)|| / ||z_0||');
    legend(CAND(:,1), 'Location', 'northeast', 'Box', 'off');
    set(gca, 'YLim', [1e-3 1]);
    save_fig(f, ['feas_' MODEL '_V'], outdir, 12, 7);

    % -- geometry : d0 and l
    f = new_fig();
    subplot(1,2,1)
    plot(d0_list*1e3, D.d0_t10/unit, 'o-', 'LineWidth', 1.5, 'Color', C.main); hold on
    vline(geo0.d0*1e3, C.grey);
    grid on; xlabel('d_0  [mm]'); ylabel(['t to 10 %  [' ulab ']']);
    subplot(1,2,2)
    plot(ls_list, D.ls_t10/unit, 'o-', 'LineWidth', 1.5, 'Color', C.main); hold on
    vline(1, C.grey);
    grid on; xlabel('trailer length, scale (1 = measured)'); ylabel(['t to 10 %  [' ulab ']']);
    save_fig(f, ['feas_' MODEL '_geom'], outdir, 14, 6);

    % -- hardware limits
    f = new_fig();
    loglog(s_list, D.hw_u1/unit, 'o-', 'LineWidth', 1.5, 'Color', C.main); hold on
    loglog(s_list, D.hw_u2/unit, 's--', 'LineWidth', 1.5, 'Color', C.second);
    vline(1, C.grey);
    set(gca, 'XTick', [0.25 0.5 1 2 4 8]);
    grid on; xlabel('limit scaled by'); ylabel(['t to 10 %  [' ulab ']']);
    legend('u_1 limit scaled', 'u_2 limit scaled', 'Location', 'northeast', 'Box', 'off');
    save_fig(f, ['feas_' MODEL '_hw'], outdir, 10, 7);
    fprintf('\nfigures and feas_%s.mat written to %s\n', MODEL, outdir);
end

% ---- helpers ----

function [b1, b2, t1, t2] = cb_wrap(Tr, z, alpha, gam, omega, k, h)
    [b1, b2, t1, t2] = Tr('control_bound', z, alpha, gam, omega, k(1), k(2), h);
end

function g = ray_gain(gmax, z0, geo, epsilon, U1, U2, cap, h_fun)
%RAY_GAIN  Largest gain admissible along the whole ray z0 -> 0:
%   amplitude at every point, and eps*gamma*h <= cap (same as gamma_ray).
    r0 = norm(z0);  zh = z0(:)/r0;  om = 2*pi/epsilon;
    g = inf;  hmax = 0;
    for r = r0*logspace(-3, 0, 200)
        g    = min(g, gmax(r*zh, geo, om, U1, U2, h_fun));
        hmax = max(hmax, h_fun(r*zh));
    end
    g = min(g, cap/(epsilon*hmax));
end

function t = interp_t(T, fr, target)
%INTERP_T  time at which ||z||/||z0|| = target, per row of T (log interpolation)
    t = zeros(size(T,1), 1);
    for c = 1:size(T,1)
        t(c) = interp1(log(fr), T(c,:), log(target));
    end
end

function [g, t10, afrac] = geo_sweep(gmax, z0of, geo0, which, vals, eps0, U1, U2, cap)
%GEO_SWEEP  Quadratic-V gain and time to 10 % when ONE geometry parameter
%   moves: which = 'd0' (values in metres) or 'l' (scale factor on l1, l2).
%   x0 stays the same physical pose, so z0 changes with the geometry.
    g = zeros(size(vals));  t10 = g;  afrac = g;
    om = 2*pi/eps0;
    for i = 1:numel(vals)
        geo = geo0;
        if strcmp(which, 'd0'), geo.d0 = vals(i);
        else, geo.l1 = geo0.l1*vals(i);  geo.l2 = geo0.l2*vals(i); end
        z = z0of(geo);
        ga = gmax(z, geo, om, U1, U2, @(zz) 1);
        g(i)   = min(ga, cap/eps0);
        t10(i) = log(10) / g(i);
        afrac(i) = g(i)/ga;                       % 1 -> amplitude-limited, <1 -> cap-limited
    end
end

function vline(x, col)
%VLINE  dotted vertical reference line in the given colour
    h = xline(x, ':');
    set(h, 'Color', col);
end

function C = palette()
    C.main   = [0.00 0.45 0.74];
    C.second = [0.85 0.33 0.10];
    C.third  = [0.47 0.67 0.19];
    C.grey   = [0.45 0.45 0.45];
    C.cat    = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56; 0.93 0.69 0.13];
end

function f = new_fig()
    f = figure('Color', 'w', 'Visible', 'on');
    set(f, 'DefaultAxesFontSize', 9, 'DefaultTextFontSize', 9);
end

function save_fig(f, name, outdir, w_cm, h_cm)
%SAVE_FIG  Vector PDF of exactly w_cm x h_cm, written with print.
    set(f, 'Units', 'centimeters', 'Position', [2 2 w_cm h_cm], ...
           'PaperUnits', 'centimeters', 'PaperSize', [w_cm h_cm], ...
           'PaperPosition', [0 0 w_cm h_cm]);
    print(f, fullfile(outdir, name), '-dpdf', '-painters');
    close(f);
end
