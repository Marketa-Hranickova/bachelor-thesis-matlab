function varargout = trailer2(task, varargin)
%TRAILER2  The whole TWO-TRAILER study (car + two off-hooked trailers, r = 4,
%   5 states) in a single file. The numbers are in params.m -- this file
%   has no numbers of its own to edit.
%
%   Not to be confused with trailers2.m, which is Victoria's original
%   script; that one is independent of this file (it only needs
%   offhooked_coeffs.m and x_to_z.m, which is why those two stay separate).
%
%   trailer2                         list of tasks, then run the default one (plots;
%                                    set DEFAULT_TASK below -- also what the Run button does)
%   trailer2('selftest')             does the code still compute what it computed before?
%   trailer2('calculations')         one closed-loop simulation, approximate model, fixed gain
%   trailer2('plots')                states, controls and convergence of that simulation
%   trailer2('compare_physical')     approximate model vs exact kinematics, as numbers
%   trailer2('compare_gain')         fixed vs adaptive gain; admissible gain and time vs epsilon
%   trailer2('compare_V')            candidate Lyapunov functions: admissible gain, time, use of the limits
%   trailer2('convergence_map')      which (epsilon, gamma) converge
%   trailer2('search_k', 15)         rank the (k12, k112, k1112) triples
%   trailer2('derive_symbolic')      Lie brackets and amplitudes symbolically (Symbolic Math Toolbox)
%
%   Every function below can be called the same way from the Command Window,
%   e.g.   [b1, b2] = trailer2('control_bound', P.z0, P.alpha, P.kappa, P.theta, ...
%                              P.gamma, P.omega, P.k12, P.k112, P.k1112, P.h_fun);
%   (the file is laid out as: tasks, simulators, controller, gains,
%   vector fields and coordinates.)
%
%   Do not put this folder and "1 trailer" on the MATLAB path together:
%   both have a params.m.

    DEFAULT_TASK = 'plots';    % what a bare  trailer2  (or the Run button) does;
                               % '' = only list the tasks

    if nargin < 1 || isempty(task)
        show_menu();
        if isempty(DEFAULT_TASK), return, end
        fprintf('Running the default task ''%s'' (change DEFAULT_TASK at the top of trailer2.m)\n\n', DEFAULT_TASK);
        task = DEFAULT_TASK;
    end
    fns = function_table();
    if ~isfield(fns, task)
        error('trailer2:unknownTask', ...
              'Unknown task "%s" -- run trailer2 without arguments to see the list.', task);
    end
    f = fns.(task);
    if nargout == 0
        f(varargin{:});
    else
        [varargout{1:nargout}] = f(varargin{:});
    end
end

function fns = function_table()
%FUNCTION_TABLE  Name -> handle of every function in this file, so that
%   trailer2('name', ...) can reach them (handles to local functions can
%   only be made inside the file itself).
    fns = struct( ...
        'selftest', @selftest, ...
        'calculations', @calculations, ...
        'plots', @plots, ...
        'compare_physical', @compare_physical, ...
        'compare_gain', @compare_gain, ...
        'compare_V', @compare_V, ...
        'convergence_map', @convergence_map, ...
        'search_k', @search_k, ...
        'derive_symbolic', @derive_symbolic, ...
        'simulate_case', @simulate_case, ...
        'simulate_case_adaptive', @simulate_case_adaptive, ...
        'simulate_case_physical', @simulate_case_physical, ...
        'rhs_approx', @rhs_approx, ...
        'u_ctrl', @u_ctrl, ...
        'osc_terms', @osc_terms, ...
        'amps', @amps, ...
        'control_bound', @control_bound, ...
        'gamma_max', @gamma_max, ...
        'eta_ray', @eta_ray, ...
        'fields_approx', @fields_approx, ...
        'fields_exact', @fields_exact, ...
        'z_to_x', @z_to_x);
end

function show_menu()
    fprintf('\ntrailer2(task, ...)   -- two trailers, r = 4, numbers in params.m\n\n');
    fprintf('  trailer2(''selftest'')           does the code still compute what it did before\n');
    fprintf('  trailer2(''calculations'')       one closed-loop simulation (approximate model)\n');
    fprintf('  trailer2(''plots'')              states / controls / convergence of that simulation\n');
    fprintf('  trailer2(''compare_physical'')   approx vs exact kinematics as numbers (slow)\n');
    fprintf('  trailer2(''compare_gain'')       fixed vs adaptive gain, gamma and time vs epsilon\n');
    fprintf('  trailer2(''compare_V'')          candidate Lyapunov functions ranked under the limits\n');
    fprintf('  trailer2(''convergence_map'')    which (epsilon, gamma) converge\n');
    fprintf('  trailer2(''search_k'', 15)       rank the (k12, k112, k1112) triples\n');
    fprintf('  trailer2(''derive_symbolic'')    symbolic derivation (needs the Symbolic Math Toolbox)\n\n');
    fprintf('  Any function in this file can be called as trailer2(''name'', args...).\n\n');
end

% ---- TASKS -- what you run ----

function selftest()
%SELFTEST  Does this folder still compute what it computed before?
%
%   Run after ANY edit to the code (about 30 s). Checks
%     1. the geometry-only numbers (alpha, kappa, theta, beta) and z0,
%     2. the control bound, gamma_max and one value of u_ctrl at z0,
%     3. a 12 s closed-loop run of all three simulators,
%   against reference values, and prints PASSED or which number moved.
%
%   The references were produced on 18 Oct 2026 with params.m as it was
%   then (gamma = 2.7498e-05, eps = 3, k = (14,5,1), sqrt V with
%   delta = 0.005, x0 = (-0.2, -0.2, 45, 45, 45 deg)). If you change those
%   ON PURPOSE, this test fails on purpose: check the new numbers, then
%   update the table below.

P = params();
n_fail = 0;
A = {P.alpha, P.kappa, P.theta};
K = {P.omega, P.k12, P.k112, P.k1112};

% ---- reference values ----------------------------------------------------
ref.coeffs = [7.355682212287460e-03  2.329616385637499e-02 ...
              6.111464836543101e-02 -1.085630195308095e-04];
ref.z0     = [-2.000000000000000e-01  7.853981633974517e-01 ...
               1.999999999999998e-01 -5.143399189943022e-02  2.492339780897482e-03];
ref.bound  = [ 2.200003688803968e-01  5.335988881155652e-01];
ref.gmax   =   2.749784700008027e-05;
ref.uctrl  = [-8.377798114535664e-03  2.109068987325568e-01];   % at t = 0.37 s

% 12 s runs; the ODE solver may differ in the last digits between MATLAB
% and Octave, hence the looser tolerance for these
ref.z_approx = [-1.999210939499e-01  7.850883000634e-01  1.999172538071e-01 ...
                -5.141390261616e-02  2.485627563645e-03];
ref.u_approx = [ 2.200003688804e-01  5.247694189120e-01];
ref.z_adapt  = [-1.999210855864e-01  7.850882672199e-01  1.999172449357e-01 ...
                -5.141390049781e-02  2.485626801677e-03];
ref.u_adapt  = [ 2.200000000000e-01  5.247685880535e-01];
ref.x_exact  = [-1.999595205010e-01 -1.999301954083e-01  7.850883015238e-01 ...
                 8.395618896470e-01  4.836880264283e-01];
ref.u_exact  = [ 2.200003688804e-01  5.247694189120e-01];

% ---- 1. geometry-only numbers --------------------------------------------
n_fail = n_fail + selftest_check('alpha, kappa, theta, beta', ...
    [P.alpha P.kappa P.theta P.beta], ref.coeffs, 1e-12);
n_fail = n_fail + selftest_check('z0 = x_to_z(x0)', P.z0(:).', ref.z0, 1e-12);

% ---- 2. controller at z0 -------------------------------------------------
[b1, b2] = control_bound(P.z0, A{:}, P.gamma, K{:}, P.h_fun);
n_fail = n_fail + selftest_check('control_bound at z0', [b1 b2], ref.bound, 1e-12);
g = gamma_max(P.z0, A{:}, K{:}, P.u1_max, P.u2_max, P.h_fun);
n_fail = n_fail + selftest_check('gamma_max at z0', g, ref.gmax, 1e-10);
[u1, u2] = u_ctrl(0.37, P.z0, K{:}, A{:}, P.gamma, P.h_fun);
n_fail = n_fail + selftest_check('u_ctrl(0.37 s) at z0', [u1 u2], ref.uctrl, 1e-12);

% ---- 3. 12 s closed loop, all three simulators ----------------------------
T = {0, 12, P.dt};                         % t_start, t_end, dt
fprintf('\nsimulating 12 s (approx / adaptive / exact) ...\n');

[~, z, u1, u2] = simulate_case(P.z0, A{:}, P.gamma, P.epsilon, K{:}, T{:}, P.h_fun);
n_fail = n_fail + selftest_check('approx   z(12 s)', z(end,:), ref.z_approx, 1e-7);
n_fail = n_fail + selftest_check('approx   max|u|', [max(abs(u1)) max(abs(u2))], ref.u_approx, 1e-7);

[~, z, u1, u2] = simulate_case_adaptive(P.z0, A{:}, P.epsilon, K{:}, T{:}, P.h_fun, ...
                                        [], [], P.u1_max, P.u2_max, P.cap_product);
n_fail = n_fail + selftest_check('adaptive z(12 s)', z(end,:), ref.z_adapt, 1e-7);
n_fail = n_fail + selftest_check('adaptive max|u|', [max(abs(u1)) max(abs(u2))], ref.u_adapt, 1e-7);

[~, x, ~, u1, u2] = simulate_case_physical(P.x0, P.d0, P.d1, P.l1, P.l2, A{:}, ...
                                           P.gamma, P.epsilon, K{:}, T{:}, P.h_fun);
n_fail = n_fail + selftest_check('exact    x(12 s)', x(end,:), ref.x_exact, 1e-7);
n_fail = n_fail + selftest_check('exact    max|u|', [max(abs(u1)) max(abs(u2))], ref.u_exact, 1e-7);

% ---- verdict ---------------------------------------------------------------
fprintf('\n');
if n_fail == 0
    fprintf('=== SELFTEST PASSED ===\n');
else
    fprintf('=== SELFTEST FAILED: %d selftest_check(s) moved -- see the FAIL lines above ===\n', n_fail);
end
end

function bad = selftest_check(label, got, want, tol)
%SELFTEST_CHECK  Compare two numeric vectors; print one line; return 1 if they differ.
    err = max(abs(got(:) - want(:)));
    ok  = err <= tol;
    if ok, tag = 'ok  '; else, tag = 'FAIL'; end
    fprintf('  [%s] %-26s max |difference| = %.2e  (tol %.0e)\n', tag, label, err, tol);
    bad = double(~ok);
end

function [t_sim, z_sim, u1_sim, u2_sim, normz_sim] = calculations()
% CALCULATIONS  One closed-loop simulation with the numbers from params.m.
%
%   Runs the approximate (nilpotent) model at the fixed gain. As a command
%   it leaves t_sim, z_sim, u1_sim, u2_sim, normz_sim in the base workspace;
%   with outputs, [t,z,u1,u2,nz] = trailer2('calculations'), it returns them.
%   No Symbolic Toolbox needed: the symbolic derivation that used to be
%   PART 1-2 of this file is now the task derive_symbolic.

    P = params();
    d0 = P.d0;  d1 = P.d1;  l1 = P.l1;  l2 = P.l2;
    alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;  beta = P.beta;
    gamma = P.gamma;  epsilon = P.epsilon;  omega = P.omega;
    k12 = P.k12;  k112 = P.k112;  k1112 = P.k1112;
    z0 = P.z0;  t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;
    U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

    fprintf('geometry : d0=%.3f  d1=%.3f  l1=%.3f  l2=%.3f  [m]\n', d0, d1, l1, l2);
    fprintf('coeffs   : alpha=%.6f  kappa=%.6f  theta=%.6f\n', alpha, kappa, theta);

    % ---- One closed-loop simulation (approximate model, fixed gain) ----
    h_fun  = P.h_fun;

    fprintf('\n========================================\n')
    fprintf('  Simulation\n')
    fprintf('  z0      = [%.4f; %.4f; %.4f; %.4f; %.4f],  ||z0|| = %.4f\n', z0, norm(z0));
    fprintf('  time    : [%.1f, %.1f], dt = %.1e\n', t_start, t_end, dt);
    fprintf('  epsilon = %.2f, omega = %.4f, gamma = %.4e\n', epsilon, omega, gamma);
    fprintf('  V       : %s\n', P.lyap_name);
    fprintf('  k       = (%d, %d, %d)\n', k12, k112, k1112);
    fprintf('========================================\n')

    [t_sim, z_sim, u1_sim, u2_sim] = simulate_case( ...
        z0, alpha, kappa, theta, gamma, ...
        epsilon, omega, k12, k112, k1112, ...
        t_start, t_end, dt, h_fun);

    fprintf('  %d steps\n', length(t_sim));
    fprintf('  final state : [%.6f; %.6f; %.6f; %.6f; %.6f]\n', z_sim(end,:));
    fprintf('  final norm  : %.6f\n', norm(z_sim(end,:)));
    fprintf('  max |u1|    : %.3f\n', max(abs(u1_sim)));
    fprintf('  max |u2|    : %.3f\n', max(abs(u2_sim)));

    normz_sim = sqrt(sum(z_sim.^2, 2));
    fprintf('  ||z(0)|| = %.6f  ->  ||z(end)|| = %.6f\n', normz_sim(1), normz_sim(end));

    fprintf('\n=== Done. ===\n')

    if nargout == 0
        assignin('base', 't_sim', t_sim);       assignin('base', 'z_sim', z_sim);
        assignin('base', 'u1_sim', u1_sim);     assignin('base', 'u2_sim', u2_sim);
        assignin('base', 'normz_sim', normz_sim);
    end
end

function plots()
% PLOTS  Figures for the single simulated case.
%
%   Runs calculations (which reads params.m) and draws three figures:
%   states, controls against the TurtleBot3 limits, convergence of ||z||.
%   Approximate model only -- for the exact kinematics see compare_physical.

    [t_sim, z_sim, u1_sim, u2_sim, normz_sim] = calculations();
    P = params();
    d0 = P.d0;  d1 = P.d1;  l1 = P.l1;  l2 = P.l2;
    U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

    % -------------------------------------------------------------- style
    col_main  = [0.00 0.20 0.70];
    col_alt   = [0.85 0.33 0.10];
    col_limit = [0.45 0.45 0.45];

    state_labels = {'z_1', 'z_2', 'z_3', 'z_4', 'z_5'};
    geom_str = sprintf('d_0=%.3f, l_1=%.3f, d_1=%.3f, l_2=%.3f m', d0, l1, d1, l2);

    % ---- FIGURE 1 -- state trajectories ----
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

    % ---- FIGURE 2 -- control inputs, against the hardware limits ----
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

    % ---- FIGURE 3 -- convergence ----
    figure('Name', 'Convergence', 'NumberTitle', 'off');
    semilogy(t_sim, max(normz_sim, realmin), 'Color', col_main, 'LineWidth', 1.5);
    grid on
    xlabel('Time [s]');
    ylabel('||z(t)||');
    title('Convergence of the state norm');

    % ---- How far outside the robot's limits is this? ----
    r1 = max(abs(u1_sim))/U1_MAX;
    r2 = max(abs(u2_sim))/U2_MAX;
    fprintf('\n========================================\n')
    fprintf('  Control demand vs TurtleBot3 Burger\n')
    fprintf('========================================\n')
    fprintf('  max |u1| = %8.4f m/s    limit %.2f    -> %6.1f %% of the limit\n', ...
            max(abs(u1_sim)), U1_MAX, 100*r1);
    fprintf('  max |u2| = %8.4f rad/s  limit %.2f    -> %6.1f %% of the limit\n', ...
            max(abs(u2_sim)), U2_MAX, 100*r2);
    if r1 > 1 + 1e-4 || r2 > 1 + 1e-4      % 1e-4: P.gamma is rounded to 5 digits
        fprintf('  -> the commanded controls do NOT fit on this robot.\n');
    else
        fprintf('  -> the commanded controls fit within the limits.\n');
    end
end

function compare_physical()
% COMPARE_PHYSICAL
% Same controller, two plants:
%   (a) the nilpotent approximation  z' = u1*g1(z) + u2*g2(z)
%       -- what simulate_case / calculations / plots show
%   (b) the exact kinematics         x' = u1*f1(x) + u2*f2(x)
%       -- what the real robot is
% The gap between them is the remainder the approximation drops, so this
% is the honest preview of the lab run.

    P = params();
    d0 = P.d0;  d1 = P.d1;  l1 = P.l1;  l2 = P.l2;
    alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;
    gamma = P.gamma;  epsilon = P.epsilon;  omega = P.omega;
    k12 = P.k12;  k112 = P.k112;  k1112 = P.k1112;
    x0 = P.x0;  z0 = P.z0;
    t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;

    h_fun  = P.h_fun;

    % ---- PART 0 -- self-check of the coordinate transform ----
    % x_to_z.m hard-codes the closed form of the frame F(0). Two checks: the
    % determinant must match eq. (4) of the paper, and z_to_x must invert it.
    fprintf('=== coordinate transform self-check ===\n');
    [~, ~, ~, ~, M0] = offhooked_coeffs(d0, d1, l1, l2);
    det_paper = -(d0 + l1)*(d0 + l2)*(d1 + l2) / (l1^4*l2^4);
    err_det   = abs(det(M0) - det_paper) / abs(det_paper);

    rng(0);
    err_rt = 0;
    for i = 1:3
        xt = [0.5*randn; 0.5*randn; 0.8*randn; 0.8*randn; 0.8*randn];
        err_rt = max(err_rt, max(abs(z_to_x(x_to_z(xt, d0, d1, l1, l2), ...
                                            d0, d1, l1, l2) - xt)));
    end
    fprintf('  rel. error in det(M0) vs eq. (4)  : %.3e\n', err_det);
    fprintf('  max |z_to_x(x_to_z(x)) - x|       : %.3e\n', err_rt);
    assert(err_det < 1e-10, 'M0 does not match the determinant in the paper');
    assert(err_rt  < 1e-9,  'z_to_x is not the inverse of x_to_z');
    disp('   PASSED')

    % ---- PART 1 -- both simulations, identical controller ----
    fprintf('\n========================================\n');
    fprintf('  x0      = [%.3f; %.3f; %.1f deg; %.1f deg; %.1f deg]\n', ...
            x0(1), x0(2), rad2deg(x0(3)), rad2deg(x0(4)), rad2deg(x0(5)));
    fprintf('  z0      = [%.4f; %.4f; %.4f; %.4f; %.4f]\n', z0);
    fprintf('  time    : [%.1f, %.1f],  epsilon = %.2f,  gamma = %.4e\n', ...
            t_start, t_end, epsilon, gamma);
    fprintf('  V       : %s\n', P.lyap_name);
    fprintf('  k       = (%d, %d, %d)\n', k12, k112, k1112);
    fprintf('========================================\n');

    fprintf('\nrunning approximate (z) model ...\n');
    [t_a, z_a, u1_a, u2_a] = simulate_case( ...
        z0, alpha, kappa, theta, gamma, ...
        epsilon, omega, k12, k112, k1112, ...
        t_start, t_end, dt, h_fun);

    fprintf('running exact (x) model ...\n');
    [t_p, x_p, z_p, u1_p, u2_p] = simulate_case_physical( ...
        x0, d0, d1, l1, l2, alpha, kappa, theta, gamma, ...
        epsilon, omega, k12, k112, k1112, ...
        t_start, t_end, dt, h_fun);

    % ---- PART 2 -- comparison ----
    fprintf('\n--- final state in z ---\n');
    fprintf('  approximate : [%9.5f %9.5f %9.5f %9.5f %9.5f]  ||z|| = %.5f\n', ...
            z_a(end,:), norm(z_a(end,:)));
    fprintf('  exact       : [%9.5f %9.5f %9.5f %9.5f %9.5f]  ||z|| = %.5f\n', ...
            z_p(end,:), norm(z_p(end,:)));
    fprintf('  ||z0||      = %.5f\n', norm(z0));

    fprintf('\n--- final state in x (what the robot physically does) ---\n');
    fprintf('  start       : pos (%.4f, %.4f) m,  phi %.2f deg, th1 %.2f deg, th2 %.2f deg\n', ...
            x0(1), x0(2), rad2deg(x0(3)), rad2deg(x0(4)), rad2deg(x0(5)));
    fprintf('  end         : pos (%.4f, %.4f) m,  phi %.2f deg, th1 %.2f deg, th2 %.2f deg\n', ...
            x_p(end,1), x_p(end,2), rad2deg(x_p(end,3)), ...
            rad2deg(x_p(end,4)), rad2deg(x_p(end,5)));
    fprintf('  ||x0|| = %.5f   ->   ||x_end|| = %.5f\n', norm(x0), norm(x_p(end,:)));

    t_c  = linspace(t_start, min(t_a(end), t_p(end)), 2000).';
    za_c = interp1(t_a, z_a, t_c);
    zp_c = interp1(t_p, z_p, t_c);
    drift = sqrt(sum((za_c - zp_c).^2, 2));
    fprintf('\n--- how far the exact model drifts from the approximation ---\n');
    fprintf('  max ||z_exact - z_approx|| : %.5f  (at t = %.2f s)\n', ...
            max(drift), t_c(find(drift == max(drift), 1)));
    fprintf('  same, relative to ||z0||   : %.1f %%\n', 100*max(drift)/norm(z0));

    fprintf('\n--- controls (unchanged by the plant, same formula both runs) ---\n');
    fprintf('  approximate : max |u1| = %.4f,  max |u2| = %.4f\n', ...
            max(abs(u1_a)), max(abs(u2_a)));
    fprintf('  exact       : max |u1| = %.4f,  max |u2| = %.4f\n', ...
            max(abs(u1_p)), max(abs(u2_p)));
    fprintf('  TurtleBot3  : u1_max  = %.4f,  u2_max  = %.4f\n', P.u1_max, P.u2_max);

    % ---- PART 3 -- plots ----
    figure('Name', 'approximate vs exact model', 'Color', 'w');

    subplot(2,2,1);
    plot(t_a, sqrt(sum(z_a.^2, 2)), 'LineWidth', 1.1); hold on
    plot(t_p, sqrt(sum(z_p.^2, 2)), 'LineWidth', 1.1);
    grid on; xlabel('t [s]'); ylabel('||z||');
    legend('approximate (z model)', 'exact (x model)', 'Location', 'best');
    title('norm in privileged coordinates');

    subplot(2,2,2);
    plot(x_p(:,1), x_p(:,2), 'LineWidth', 1.1); hold on
    plot(x_p(1,1), x_p(1,2), 'o', 'MarkerFaceColor', 'auto');
    plot(0, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);
    grid on; axis equal; xlabel('x_1 [m]'); ylabel('x_2 [m]');
    legend('path', 'start', 'goal', 'Location', 'best');
    title('where the vehicle actually drives');

    subplot(2,2,3);
    plot(t_p, rad2deg(x_p(:,3:5)), 'LineWidth', 1.1);
    grid on; xlabel('t [s]'); ylabel('angle [deg]');
    legend('\phi', '\theta_1', '\theta_2', 'Location', 'best');
    title('physical angles, exact model');

    subplot(2,2,4);
    plot(t_c, drift, 'LineWidth', 1.1);
    grid on; xlabel('t [s]'); ylabel('||z_{exact} - z_{approx}||');
    title('cost of the nilpotent approximation');

    fprintf('\n=== Done. ===\n');
end

function compare_gain()
% COMPARE_GAIN  Fixed vs. adaptive gain, and how the admissible gain
% depends on epsilon.
%
%   fixed    : gamma = P.gamma for the whole run.
%   adaptive : gamma recomputed every interval = largest value that still
%              respects the robot's |u1|,|u2| limits (gamma_max),
%              additionally capped so eps*gamma*h(z) <= P.cap_product

    clc; close all;

    P = params();
    alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;
    k12 = P.k12;  k112 = P.k112;  k1112 = P.k1112;
    z0 = P.z0;  t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;
    U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

    h_fun  = P.h_fun;

    fprintf('geometry : d0=%.3f  d1=%.3f  l1=%.3f  l2=%.3f  [m]\n', P.d0, P.d1, P.l1, P.l2);
    fprintf('z0 = [%g %g %g %g %g],  eps = %g,  k = (%d,%d,%d)\n', z0, P.epsilon, k12, k112, k1112);
    fprintf('V  : %s\n', P.lyap_name);

    % ---------------- fixed gain ------------------------------------------------
    fprintf('\n--- fixed gain, gamma = %g ---\n', P.gamma);
    [t_f, z_f, u1_f, u2_f] = simulate_case(z0, alpha, kappa, theta, P.gamma, ...
        P.epsilon, P.omega, k12, k112, k1112, t_start, t_end, dt, h_fun);
    fprintf('  final norm : %.6f\n', norm(z_f(end,:)));
    fprintf('  max |u1|   : %.3f  (limit %.2f)\n', max(abs(u1_f)), U1_MAX);
    fprintf('  max |u2|   : %.3f  (limit %.2f)\n', max(abs(u2_f)), U2_MAX);

    % ---------------- adaptive gain ---------------------------------------------
    fprintf('\n--- adaptive gain, cap eps*gamma*h <= %g ---\n', P.cap_product);
    [t_a, z_a, u1_a, u2_a, gamma_a] = simulate_case_adaptive(z0, alpha, kappa, theta, ...
        P.epsilon, P.omega, k12, k112, k1112, t_start, t_end, dt, h_fun, [], [], ...
        U1_MAX, U2_MAX, P.cap_product);
    fprintf('  final norm : %.6f\n', norm(z_a(end,:)));
    fprintf('  max |u1|   : %.3f  (limit %.2f)\n', max(abs(u1_a)), U1_MAX);
    fprintf('  max |u2|   : %.3f  (limit %.2f)\n', max(abs(u2_a)), U2_MAX);
    fprintf('  gamma used : min %.4g   max %.4g   (fixed run used %.4g throughout)\n', ...
            min(gamma_a), max(gamma_a), P.gamma);

    % ---------------- how long would the adaptive controller actually take? ----
    gamma_bar = mean(gamma_a);
    targets   = [0.5 0.1 0.01 0.001];
    t_need    = eta_ray(z0, gamma_bar, h_fun, targets);

    fprintf('\n--- estimated time to converge at the mean adaptive gain ---\n');
    fprintf('  mean gamma used : %.4g   (gradient-flow estimate, eta_ray)\n', gamma_bar);
    for i = 1:numel(targets)
        fprintf('  ||z|| down to %5.1f%% of ||z0|| : %10.0f s = %7.2f h\n', ...
                targets(i)*100, t_need(i), t_need(i)/3600);
    end
    fprintf('  (the run above lasted %.0f s, so the adaptive gain only changes\n', t_end - t_start);
    fprintf('   this estimate if it actually moved away from the fixed gamma)\n');

    col_f = [0.00 0.20 0.70];
    col_a = [0.85 0.33 0.10];

    eps_sweep    = [0.05 0.1 0.2 0.5 1 2 4 8 16 32 64 128];
    n_sweep      = numel(eps_sweep);
    gmax_sweep   = zeros(1, n_sweep);
    hrs10_sweep  = nan(1, n_sweep);
    hrs1_sweep   = nan(1, n_sweep);

    % largest h on the ray (for the contraction check eps*gamma*h <= cap)
    r0   = norm(z0);
    hmax = max(arrayfun(@(r) h_fun(r*z0/r0), r0*logspace(-3, 0, 200)));

    fprintf('\n===== gamma_max and time-to-converge vs epsilon =====\n');
    fprintf('(gamma_max at z0 for the velocity limits; time from eta_ray;\n');
    fprintf(' "--" = eps*gamma*h exceeds cap_product = %.2f, estimate not valid)\n', P.cap_product);
    fprintf('%10s %14s %14s %14s %14s\n', 'eps [s]', 'gamma_max', 'eps*g*hmax', 'hours to 10%', 'hours to 1%');
    for i = 1:n_sweep
        Q = params(eps_sweep(i));
        gmax_sweep(i) = gamma_max(z0, alpha, kappa, theta, Q.omega, ...
                                  k12, k112, k1112, U1_MAX, U2_MAX, h_fun);

        egh = Q.epsilon * gmax_sweep(i) * hmax;
        if egh <= P.cap_product
            t = eta_ray(z0, gmax_sweep(i), h_fun, [0.1 0.01]) / 3600;
            hrs10_sweep(i) = t(1);
            hrs1_sweep(i)  = t(2);
            fprintf('%10.2f %14.3e %14.3e %14.2f %14.2f\n', ...
                    Q.epsilon, gmax_sweep(i), egh, hrs10_sweep(i), hrs1_sweep(i));
        else
            fprintf('%10.2f %14.3e %14.3e %14s %14s\n', ...
                    Q.epsilon, gmax_sweep(i), egh, '--', '--');
        end
    end
    fprintf(['NOTE: large eps looks fast here but the nilpotent approximation stops\n' ...
             '      describing the real kinematics beyond eps ~ 4 s (see params.m) --\n' ...
             '      check any retune with compare_physical.\n']);

    figure('Name', 'gamma_max and ETA vs epsilon', 'NumberTitle', 'off');
    subplot(1,2,1)
    loglog(eps_sweep, gmax_sweep, 'o-', 'Color', col_a, 'LineWidth', 1.5);
    grid on; xlabel('\epsilon  [s]'); ylabel('\gamma_{max}  (hardware-safe)');
    title('Largest hardware-safe gain vs. period');

    subplot(1,2,2)
    loglog(eps_sweep, hrs10_sweep, 'o-',  'Color', col_a, 'LineWidth', 1.5); hold on
    loglog(eps_sweep, hrs1_sweep,  's--', 'Color', col_a, 'LineWidth', 1.0);
    grid on; xlabel('\epsilon  [s]'); ylabel('hours to converge');
    legend('to 10% of ||z_0||', 'to 1% of ||z_0||', 'Location', 'best');
    title('Estimated time to converge at gamma_{max}(\epsilon)');

    % ---------------- plot: ||z(t)|| --------------------------------------------
    figure('Name', 'Fixed vs adaptive gain -- convergence', 'NumberTitle', 'off');
    semilogy(t_f, vecnorm(z_f, 2, 2), 'Color', col_f, 'LineWidth', 1.5); hold on
    semilogy(t_a, vecnorm(z_a, 2, 2), 'Color', col_a, 'LineWidth', 1.5);
    grid on; xlabel('t [s]'); ylabel('||z(t)||');
    legend('fixed gain', 'adaptive gain', 'Location', 'best');
    title('Convergence: fixed vs. adaptive gain');

    % ---------------- plot: controls, both cases --------------------------------
    figure('Name', 'Fixed vs adaptive gain -- controls', 'NumberTitle', 'off');
    subplot(2,1,1)
    plot(t_f, u1_f, 'Color', col_f, 'LineWidth', 1.0); hold on
    plot(t_a, u1_a, 'Color', col_a, 'LineWidth', 1.0);
    yline(U1_MAX, 'k--'); yline(-U1_MAX, 'k--');
    grid on; ylabel('u_1 [m/s]');
    legend('fixed', 'adaptive', 'limit', 'Location', 'best');
    title('Control inputs: fixed vs. adaptive gain');

    subplot(2,1,2)
    plot(t_f, u2_f, 'Color', col_f, 'LineWidth', 1.0); hold on
    plot(t_a, u2_a, 'Color', col_a, 'LineWidth', 1.0);
    yline(U2_MAX, 'k--'); yline(-U2_MAX, 'k--');
    grid on; xlabel('t [s]'); ylabel('u_2 [rad/s]');
    legend('fixed', 'adaptive', 'limit', 'Location', 'best');

    % ---------------- plot: gamma(t) for the adaptive run -----------------------
    figure('Name', 'Adaptive gain over time', 'NumberTitle', 'off');
    plot(t_a, gamma_a, 'Color', col_a, 'LineWidth', 1.5); hold on
    yline(P.gamma, 'b--');
    grid on; xlabel('t [s]'); ylabel('\gamma used');
    legend('adaptive \gamma(t)', 'fixed \gamma', 'Location', 'best');
    title('Gain chosen every interval (hardware-limited, capped \epsilon\gamma h <= cap)');

    fprintf('\n=== Done. ===\n');
end

function compare_V()
% COMPARE_V  Candidate Lyapunov functions, ranked under the robot's limits.
%
% For a radially symmetric V = phi(r), r = ||z||, the whole function enters
% the controller through one scalar factor (see params.m):
%
%       a(z) = -gamma*h(z)*G^-1(z)*z,      h(z) = phi'(r)/(gamma*r)
%
% so the candidates below differ ONLY in h. Two consequences used here:
%
%   1. Every radial V drives the state along the SAME straight ray to the
%      origin in z; only the speed along it changes. So the admissible gain
%      can be found by walking that ray.
%   2. The closed loop approximates the gradient flow  dr/dt = -gamma*h(r)*r,
%      a 1-D ODE, so the convergence time follows without simulating hours
%      of oscillation (eta_ray).
%
% Per candidate, gamma is the largest value that keeps
%   (1) the control bound within the limits at EVERY point of the ray, and
%   (2) eps*gamma*h <= P.cap_product on the ray (contraction, see params.m).
% The row marked "in use" is whatever params.m currently sets.

    P = params();
    alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;
    omega = P.omega;  k12 = P.k12;  k112 = P.k112;  k1112 = P.k1112;
    U1 = P.u1_max;  U2 = P.u2_max;
    z0 = P.z0;  r0 = norm(z0);  zhat = z0/r0;

    CAND = { ...
        ['in use: ' P.lyap_name],               P.h_fun; ...
        'quadratic  V = gamma*r^2/2',           @(z) 1; ...
        'log c=1    V = gamma/(2c)log(1+cr^2)', @(z) 1/(1 +   1*norm(z)^2); ...
        'log c=10   (Grushkovskaya''s start)',  @(z) 1/(1 +  10*norm(z)^2); ...
        'log c=100',                            @(z) 1/(1 + 100*norm(z)^2); ...
        'sqrt d=0.01 V = gamma(sqrt(r^2+d^2)-d)', @(z) 1/sqrt(norm(z)^2 + 0.01^2); ...
        'sqrt d=0.005',                         @(z) 1/sqrt(norm(z)^2 + 0.005^2) };

    rs = r0 * logspace(-3, 0, 200);

    fprintf('||z0|| = %.5f,  limits |u1|<=%.2f, |u2|<=%.2f,  k=(%d,%d,%d), eps=%.1f, cap=%.2f\n\n', ...
            r0, U1, U2, k12, k112, k1112, P.epsilon, P.cap_product);
    nC  = size(CAND, 1);
    res = zeros(nC, 7);     % gamma, t10, t1, t01, eps*g*hmax, max u1 use, max u2 use
    for c = 1:nC
        h = CAND{c, 2};

        % (1) largest gain admissible at EVERY point of the ray
        g = inf;
        hmax = 0;
        for r = rs
            g    = min(g, gamma_max(r*zhat, alpha, kappa, theta, omega, ...
                                    k12, k112, k1112, U1, U2, h));
            hmax = max(hmax, h(r*zhat));
        end
        % (2) contraction cap
        g = min(g, P.cap_product/(P.epsilon*hmax));

        % how much of each limit the chosen gain uses (worst point of the ray)
        use1 = 0;  use2 = 0;
        for r = rs
            [b1, b2] = control_bound(r*zhat, alpha, kappa, theta, g, omega, ...
                                     k12, k112, k1112, h);
            use1 = max(use1, b1/U1);  use2 = max(use2, b2/U2);
        end

        t = eta_ray(z0, g, h, [0.10 0.01 0.001]) / 3600;      % hours
        res(c, :) = [g, t, P.epsilon*g*hmax, use1, use2];
    end
    base = res(2, 3);                    % row 2 = quadratic V, time to 1%

    fprintf('%-40s %10s %8s %8s %8s %8s %9s %6s %6s\n', 'Lyapunov function', 'gamma', ...
            't->10%', 't->1%', 't->0.1%', 'speedup', 'eps*g*h', 'u1 %', 'u2 %');
    for c = 1:nC
        fprintf('%-40s %10.3e %6.2f h %6.2f h %6.2f h %7.2fx %9.2e %5.0f%% %5.0f%%\n', ...
                CAND{c,1}, res(c,1), res(c,2), res(c,3), res(c,4), base/res(c,3), ...
                res(c,5), 100*res(c,6), 100*res(c,7));
    end
    fprintf(['\nspeedup = time to 1%% of the quadratic V divided by this row''s.\n' ...
             'eps*g*h = largest eps*gamma*h on the ray (contraction, must stay <= cap).\n' ...
             'u1 %%, u2 %% = largest share of each velocity limit used on the ray;\n' ...
             '100%% marks the limit that sets gamma.\n']);
end

function convergence_map()
% CONVERGENCE_MAP  For which (epsilon, gamma) does the controller converge?
%
%   Sweeps a grid of sampling periods epsilon and gains gamma, runs a short
%   simulation of the APPROXIMATE model for each pair and classifies the
%   outcome: converges / marginal / diverges. The result is compared with
%   the proof's linear bound eps*gamma = 1 (exact hit) and 2 (breakdown);
%   that comparison is where P.cap_product = 0.8 in params.m comes from.
%
%   Uses the Lyapunov function from params.m. Takes a few minutes
%   (144 short simulations); shrink eps_list / gamma_list for a quick look.

    clc; close all;

    P = params();
    alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;
    k12 = P.k12;  k112 = P.k112;  k1112 = P.k1112;
    z0 = P.z0;  t_start = P.t_start;

    h_fun  = P.h_fun;

    % ---------------- sweep settings (fast, short horizon) --------------------
    t_end_map = 8;
    dt_map    = 1e-3;

    eps_list   = logspace(log10(0.05), log10(16), 12);
    gamma_list = logspace(log10(0.01), log10(3),  12);

    n0 = norm(z0);
    % verdict: 1 = converges, 2 = marginal, 3 = diverges / blows up
    verdict = zeros(numel(eps_list), numel(gamma_list));

    fprintf('Sweeping %d x %d = %d combinations, t_end = %g s ...\n', ...
            numel(eps_list), numel(gamma_list), numel(eps_list)*numel(gamma_list), t_end_map);
    tic
    for i = 1:numel(eps_list)
        Q       = params(eps_list(i));
        eps_i   = Q.epsilon;
        omega_i = Q.omega;
        for j = 1:numel(gamma_list)
            g_j = gamma_list(j);
            try
                [~, z_sim] = simulate_case(z0, alpha, kappa, theta, g_j, ...
                    eps_i, omega_i, k12, k112, k1112, ...
                    t_start, t_end_map, dt_map, h_fun);
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

    % ---------------- numeric check: real boundary vs. eps*gamma = 1 and 2 -----
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
end

function T = search_k(kmax, eps_override)
%SEARCH_K  All C1-admissible multiplier triples, ranked by admissible gain.
%
%   T = search_k(15)        search k = 1..15, epsilon from params.m
%   T = search_k(15, 8.0)   with epsilon = 8 s
%
%   For every triple that satisfies condition C1 it computes gamma_max (the
%   largest gain within the velocity limits at P.z0) and max_u1 (worst-case
%   |u1| at the CURRENT P.gamma). Returns the table sorted by gamma_max and
%   prints the best ten next to the paper's (8,6,5). What usually pays is a SMALL
%   k1112: the deepest bracket term carries omega^(3/4) and a factor 6
%   in u2, so it dominates the demand.

    if nargin < 1, kmax = 30; end
    P = params();
    if nargin >= 2
        P = params(eps_override);
    end
    alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;
    gamma = P.gamma;  epsilon = P.epsilon;  w = P.omega;
    z0 = P.z0;  U1_MAX = P.u1_max;  U2_MAX = P.u2_max;
    h_fun  = P.h_fun;
    % (the loop below supplies its own multipliers k)

    fprintf('geometry : d0=%.3f  d1=%.3f  l1=%.3f  l2=%.3f  [m]\n', P.d0, P.d1, P.l1, P.l2);
    fprintf('coeffs   : alpha=%.6f  kappa=%.6f  theta=%.6f\n', alpha, kappa, theta);
    fprintf('V        : %s,  max_u1 column at gamma = %.4e\n', P.lyap_name, gamma);

    rows = zeros(0, 5);
    for a = 1:kmax
        for b = 1:kmax
            for c = 1:kmax
                if satisfies_c1(a, b, c)
                    g  = gamma_max(z0, alpha, kappa, theta, w, a, b, c, ...
                                   U1_MAX, U2_MAX, h_fun);
                    u1 = control_bound(z0, alpha, kappa, theta, gamma, w, ...
                                       a, b, c, h_fun);
                    rows(end+1, :) = [a, b, c, g, u1];   %#ok<AGROW>
                end
            end
        end
    end

    [~, idx] = sort(rows(:, 4), 'descend');
    rows = rows(idx, :);

    T = array2table(rows, 'VariableNames', ...
        {'k12', 'k112', 'k1112', 'gamma_max', 'max_u1'});

    base = gamma_max(z0, alpha, kappa, theta, w, 8, 6, 5, ...
                     U1_MAX, U2_MAX, h_fun);
    fprintf('\n %d admissible triples with k <= %d,  eps = %.2f s\n', ...
            size(rows, 1), kmax, epsilon);
    fprintf(' the paper''s (8,6,5) gives gamma_max = %.3e\n\n', base);
    fprintf('%6s %6s %8s %14s %12s %10s\n', ...
            'k12', 'k112', 'k1112', 'gamma_max', 'max|u1|', 'vs paper');
    for i = 1:min(10, size(rows, 1))
        fprintf('%6d %6d %8d %14.4e %12.3f %9.2fx\n', ...
                rows(i,1), rows(i,2), rows(i,3), rows(i,4), rows(i,5), ...
                rows(i,4)/base);
    end
    fprintf('\n');
end

function ok = satisfies_c1(k12, k112, k1112)
%SATISFIES_C1  Condition C1 of Theorem 1 (Zuyev & Grushkovskaya, 2403.03341).
%
%   C1) there are no resonances of order up to 4 in each of the tuples
%       (k12, k112, k1112), (k12, 2*k112, k1112),
%       (k12, k112, 3*k1112), (k12, 2*k112, 3*k1112)
%
%   A tuple (k_1,...,k_s) is in resonance of order K if there exist SETWISE
%   COPRIME integers (c_1,...,c_s) with sum|c_i| = K and sum c_i*k_i = 0.
%
%   Smaller multipliers mean smaller commanded velocities, so it is tempting
%   to shrink them -- this is the condition that says how far you may go.

    tuples = [k12, k112,   k1112;
              k12, 2*k112, k1112;
              k12, k112,   3*k1112;
              k12, 2*k112, 3*k1112];

    ok = true;
    for r = 1:size(tuples, 1)
        if has_resonance(tuples(r, :), 4)
            ok = false;
            return
        end
    end
end

function res = has_resonance(kt, max_order)
%HAS_RESONANCE  True if the tuple has a resonance of order <= max_order.
    res = false;
    for c1 = -max_order:max_order
        for c2 = -max_order:max_order
            for c3 = -max_order:max_order
                K = abs(c1) + abs(c2) + abs(c3);
                if K == 0 || K > max_order
                    continue
                end
                if gcd(gcd(abs(c1), abs(c2)), abs(c3)) ~= 1
                    continue        % not setwise coprime
                end
                if c1*kt(1) + c2*kt(2) + c3*kt(3) == 0
                    res = true;
                    return
                end
            end
        end
    end
end

function derive_symbolic()
% DERIVE_SYMBOLIC  Lie brackets, controllability and the amplitudes a(z).
%
%   Needs the Symbolic Math Toolbox and takes a while (simplify). Run it to
%   see or re-check the algebra behind the controller; simulations do not
%   need it. (This was PART 1 and 2 of calculations.)
%
%   Part 1  vector fields g1, g2 of the nilpotent approximation, their Lie
%           brackets up to depth 4, the matrix
%               G(z) = [g1, g2, [g1,g2], [g1,[g1,g2]], [g1,[g1,[g1,g2]]]]
%           and det G (= 1, so G is invertible everywhere).
%   Part 2  a(z) = -G^-1(z)*grad Q(z) for the quadratic Q = gamma/2*||z||^2,
%           checked against amps. Any radial V gives the same formulas
%           times the scalar h(z) (see params.m).

    P = params();
    alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;  gamma = P.gamma;

    fprintf('geometry : d0=%.3f  d1=%.3f  l1=%.3f  l2=%.3f  [m]\n', P.d0, P.d1, P.l1, P.l2);
    fprintf('coeffs   : alpha=%.6f  kappa=%.6f  theta=%.6f\n', alpha, kappa, theta);

    % ---- PART 1 -- symbolic system, Lie brackets, controllability ----
    syms z1 z2 z3 z4 z5 real
    syms alpha_s kappa_s theta_s real
    syms gamma_s positive

    z_sym = [z1; z2; z3; z4; z5];

    % Vector fields of the nilpotent approximation. The system is
    %     z' = u1*g1(z) + u2*g2(z)
    % (numerically the same fields as fields_approx)
    g1_sym = [1; 0; -z2; -z3; kappa_s*z2^3 - z4];
    g2_sym = [0; 1; 0; -alpha_s*z2^2; -theta_s*z2*z3];

    fprintf('\n=== Vector field g1(z) ===\n'); disp(g1_sym)
    fprintf('\n=== Vector field g2(z) ===\n'); disp(g2_sym)

    dg1dz = jacobian(g1_sym, z_sym);
    dg2dz = jacobian(g2_sym, z_sym);

    % [g1, g2]
    Lie_g1g2_sym = simplify(dg2dz*g1_sym - dg1dz*g2_sym);
    fprintf('\n=== [g1, g2] ===\n'); disp(Lie_g1g2_sym)

    % [g1, [g1, g2]]
    dLie_g1g2_dz    = jacobian(Lie_g1g2_sym, z_sym);
    Lie_g1_g1g2_sym = simplify(dLie_g1g2_dz*g1_sym - dg1dz*Lie_g1g2_sym);
    fprintf('\n=== [g1, [g1, g2]] ===\n'); disp(Lie_g1_g1g2_sym)

    % [g1, [g1, [g1, g2]]]
    dLie_g1_g1g2_dz    = jacobian(Lie_g1_g1g2_sym, z_sym);
    Lie_g1_g1_g1g2_sym = simplify(dLie_g1_g1g2_dz*g1_sym - dg1dz*Lie_g1_g1g2_sym);
    fprintf('\n=== [g1, [g1, [g1, g2]]] ===\n'); disp(Lie_g1_g1_g1g2_sym)

    % The five brackets span the tangent space -> the system is controllable.
    G_sym = [g1_sym, g2_sym, Lie_g1g2_sym, Lie_g1_g1g2_sym, Lie_g1_g1_g1g2_sym];
    fprintf('\n=== G matrix ===\n'); disp(G_sym)

    detG = simplify(det(G_sym));
    fprintf('\n=== det(G) ===\n'); disp(detG)

    % ---- PART 2 -- the amplitudes a(z), for h = 1 ----
    % Lyapunov function Q(z) = gamma/2 * ||z||^2   (eq. 17)
    Q_sym     = gamma_s/2 * (z1^2 + z2^2 + z3^2 + z4^2 + z5^2);
    gradQ_sym = jacobian(Q_sym, z_sym).';
    fprintf('\n=== grad Q(z) ===\n'); disp(gradQ_sym)

    G_inv_sym = simplify(inv(G_sym));

    % a(z) = -h * G^{-1} * grad Q,  with h = 1
    a_sym = simplify(-G_inv_sym*gradQ_sym);

    fprintf('\n=== a(z) components, h = 1 ===\n')
    labels = {'a1   ', 'a2   ', 'a12  ', 'a112 ', 'a1112'};
    for i = 1:5
        fprintf('%s = ', labels{i}); disp(a_sym(i))
    end

    % ------------------------------------------------------------------ check
    % Every simulation uses the closed forms in amps. Check them against the
    % symbolic result at a generic point, so a typo in either place is caught.
    z_test = [0.37; -0.21; 0.44; -0.13; 0.29];

    a_from_sym = double(subs(a_sym, ...
        [z_sym.', gamma_s, alpha_s, kappa_s, theta_s], ...
        [z_test.', gamma, alpha, kappa, theta]));

    % amps with h = 1, to match the quadratic V derived above
    [b1, b2, b12, b112, b1112] = amps(z_test, alpha, kappa, theta, gamma, @(z) 1);
    a_from_amps = [b1; b2; b12; b112; b1112];

    res = max(abs(a_from_sym - a_from_amps));
    fprintf('\n=== symbolic vs amps: max |difference| = %.3e ===\n', res);
    assert(res < 1e-10, 'a(z) mismatch -- symbolic result and amps disagree');
    disp('   PASSED')

    fprintf('\n=== Done. ===\n')
end

% ---- SIMULATORS -- the three closed loops (approximate / adaptive gain / exact kinematics) ----

function [t_out, z_out, u1_out, u2_out] = simulate_case( ...
        z0, alpha, kappa, theta, gamma_val, ...
        epsilon, omega, k12, k112, k1112, ...
        t_start, t_end, dt, h_fun, g1_fun, g2_fun)
%SIMULATE_CASE  Closed loop on the NILPOTENT APPROXIMATION
%   z' = u1*g1(z) + u2*g2(z), fixed gain gamma_val.
%
%   Every epsilon seconds the state is frozen, the controller amplitudes
%   are computed from it (u_ctrl), and ode45 integrates one interval.
%
%   g1_fun, g2_fun are optional: leave them out (or pass []) and the fields
%   come from fields_approx.
%   For the exact kinematics use simulate_case_physical, for the
%   state-dependent gain simulate_case_adaptive.

    if nargin < 15, g1_fun = []; end
    if nargin < 16, g2_fun = []; end

    % the -1e-9 stops float noise (e.g. 3.0000000000000004) adding an empty interval
    N_int   = ceil((t_end - t_start)/epsilon - 1e-9);
    t_out   = [];
    z_out   = [];
    z_cur   = z0(:);
    frozen  = zeros(5, N_int);          % the frozen state of each interval
    t_edges = zeros(1, N_int);
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', dt);

    for j = 1:N_int
        t_j   = t_start + (j-1)*epsilon;
        t_jp1 = min(t_j + epsilon, t_end);
        if t_j >= t_end
            frozen = frozen(:, 1:j-1);
            t_edges = t_edges(1:j-1);
            break
        end

        z_frozen     = z_cur;
        frozen(:, j) = z_frozen;
        t_edges(j)   = t_j;

        odefun = @(t, z) rhs_approx(t, z, z_frozen, g1_fun, g2_fun, omega, ...
                               k12, k112, k1112, alpha, kappa, theta, ...
                               gamma_val, h_fun);

        % output grid of this interval; it must END exactly at t_jp1,
        % otherwise the next interval starts from a state at an earlier time
        tspan_j = t_j:dt:t_jp1;
        if tspan_j(end) < t_jp1 - 1e-12
            tspan_j = [tspan_j, t_jp1];
        end
        if numel(tspan_j) < 2
            tspan_j = [t_j, t_jp1];
        end

        [t_seg, z_seg] = ode45(odefun, tspan_j, z_cur, options);

        if isempty(t_out)
            t_out = t_seg;
            z_out = z_seg;
        else
            t_out = [t_out; t_seg(2:end)];
            z_out = [z_out; z_seg(2:end, :)];
        end
        z_cur = z_seg(end, :).';
    end

    n      = numel(t_out);
    u1_out = zeros(n, 1);
    u2_out = zeros(n, 1);
    n_int  = size(frozen, 2);
    for i = 1:n
        j_idx = min(max(floor((t_out(i) - t_start)/epsilon) + 1, 1), n_int);
        zf    = frozen(:, j_idx);
        [u1_out(i), u2_out(i)] = u_ctrl(t_out(i), zf, omega, k12, k112, ...
                                        k1112, alpha, kappa, theta, gamma_val, h_fun);
    end
end

function [t_out, z_out, u1_out, u2_out, gamma_out] = simulate_case_adaptive( ...
        z0, alpha, kappa, theta, ...
        epsilon, omega, k12, k112, k1112, ...
        t_start, t_end, dt, h_fun, g1_fun, g2_fun, ...
        u1_max, u2_max, cap_product)
%SIMULATE_CASE_ADAPTIVE  Like simulate_case, but the gain is chosen anew
%   on every sampling interval: the largest gamma that keeps |u1|, |u2|
%   within (u1_max, u2_max) at the frozen state (gamma_max), capped so
%   that eps*gamma*h(z) <= cap_product.
%
%   g1_fun, g2_fun may be [] -- then the fields come from fields_approx.
%   gamma_out is the gain in use at every output time.

    % the -1e-9 stops float noise (e.g. 3.0000000000000004) adding an empty interval
    N_int     = ceil((t_end - t_start)/epsilon - 1e-9);
    t_out     = [];
    z_out     = [];
    z_cur     = z0(:);
    frozen    = zeros(5, N_int);        % frozen state of each interval
    gamma_int = zeros(1, N_int);        % gain actually used on each interval
    options   = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', dt);

    for j = 1:N_int
        t_j   = t_start + (j-1)*epsilon;
        t_jp1 = min(t_j + epsilon, t_end);
        if t_j >= t_end
            frozen    = frozen(:, 1:j-1);
            gamma_int = gamma_int(1:j-1);
            break
        end

        z_frozen     = z_cur;
        frozen(:, j) = z_frozen;

        % ---- adaptive gain for this interval -------------------------------
        % h_fun must go in here too: the gain is chosen from the same
        % amplitudes the controller then uses, so a mismatch would pick a
        % gain that violates the limits by exactly the factor h.
        g_hw = gamma_max(z_frozen, alpha, kappa, theta, omega, ...
                          k12, k112, k1112, u1_max, u2_max, h_fun);
        % contraction cap eps*gamma*h <= cap_product (params.m, rule (2)).
        % For the quadratic V h = 1 and this is the old eps*gamma <= cap;
        % for the sqrt / log V h varies with the state and must be included.
        g_j  = min(g_hw, cap_product/(epsilon*h_fun(z_frozen)));
        gamma_int(j) = g_j;

        odefun = @(t, z) rhs_approx(t, z, z_frozen, g1_fun, g2_fun, omega, ...
                               k12, k112, k1112, alpha, kappa, theta, ...
                               g_j, h_fun);

        % output grid of this interval; it must END exactly at t_jp1,
        % otherwise the next interval starts from a state at an earlier time
        tspan_j = t_j:dt:t_jp1;
        if tspan_j(end) < t_jp1 - 1e-12
            tspan_j = [tspan_j, t_jp1];
        end
        if numel(tspan_j) < 2
            tspan_j = [t_j, t_jp1];
        end

        [t_seg, z_seg] = ode45(odefun, tspan_j, z_cur, options);

        if isempty(t_out)
            t_out = t_seg;
            z_out = z_seg;
        else
            t_out = [t_out; t_seg(2:end)];
            z_out = [z_out; z_seg(2:end, :)];
        end
        z_cur = z_seg(end, :).';
    end

    n         = numel(t_out);
    u1_out    = zeros(n, 1);
    u2_out    = zeros(n, 1);
    gamma_out = zeros(n, 1);
    n_int     = size(frozen, 2);
    for i = 1:n
        j_idx = min(max(floor((t_out(i) - t_start)/epsilon) + 1, 1), n_int);
        zf    = frozen(:, j_idx);
        g_j   = gamma_int(j_idx);
        gamma_out(i) = g_j;
        [u1_out(i), u2_out(i)] = u_ctrl(t_out(i), zf, omega, k12, k112, ...
                                        k1112, alpha, kappa, theta, g_j, h_fun);
    end
end

function [t_out, x_out, z_out, u1_out, u2_out] = simulate_case_physical( ...
        x0, d0, d1, l1, l2, alpha, kappa, theta, gamma_val, ...
        epsilon, omega, k12, k112, k1112, ...
        t_start, t_end, dt, h_fun)
%SIMULATE_CASE_PHYSICAL  Closed loop on the REAL model.
%
%   simulate_case integrates the nilpotent approximation
%       z' = u1*g1(z) + u2*g2(z),
%   i.e. the truncated system the controller was designed on. This file
%   integrates the exact kinematics instead,
%       x' = u1*f1(x) + u2*f2(x),
%   while driving it with exactly the same controller. The difference
%   between the two runs is the remainder the approximation throws away --
%   which is what the real robot will actually feel.
%
%   Loop structure per sampling interval (also the loop the robot runs):
%       1. read the physical state x
%       2. z = x_to_z(x)                 <- state conversion
%       3. u1, u2 from the controller at that frozen z
%       4. apply u1, u2 to the plant     <- NO conversion back: u1, u2 are
%                                           already the physical driving
%                                           and steering velocities

    % the -1e-9 stops float noise (e.g. 3.0000000000000004) adding an empty interval
    N_int   = ceil((t_end - t_start)/epsilon - 1e-9);
    t_out   = [];
    x_out   = [];
    x_cur   = x0(:);
    frozen  = zeros(5, N_int);          % frozen z of each interval
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', dt);

    for j = 1:N_int
        t_j   = t_start + (j-1)*epsilon;
        t_jp1 = min(t_j + epsilon, t_end);
        if t_j >= t_end
            frozen = frozen(:, 1:j-1);
            break
        end

        % The controller speaks z, the plant speaks x.
        z_frozen     = x_to_z(x_cur, d0, d1, l1, l2);
        frozen(:, j) = z_frozen;

        odefun = @(t, x) rhs_x(t, x, z_frozen, d0, d1, l1, l2, omega, ...
                               k12, k112, k1112, alpha, kappa, theta, ...
                               gamma_val, h_fun);

        % output grid of this interval; it must END exactly at t_jp1,
        % otherwise the next interval starts from a state at an earlier time
        tspan_j = t_j:dt:t_jp1;
        if tspan_j(end) < t_jp1 - 1e-12
            tspan_j = [tspan_j, t_jp1];
        end
        if numel(tspan_j) < 2
            tspan_j = [t_j, t_jp1];
        end

        [t_seg, x_seg] = ode45(odefun, tspan_j, x_cur, options);

        if isempty(t_out)
            t_out = t_seg;
            x_out = x_seg;
        else
            t_out = [t_out; t_seg(2:end)];
            x_out = [x_out; x_seg(2:end, :)];
        end
        x_cur = x_seg(end, :).';
    end

    n      = numel(t_out);
    u1_out = zeros(n, 1);
    u2_out = zeros(n, 1);
    n_int  = size(frozen, 2);
    for i = 1:n
        j_idx = min(max(floor((t_out(i) - t_start)/epsilon) + 1, 1), n_int);
        [u1_out(i), u2_out(i)] = u_ctrl(t_out(i), frozen(:, j_idx), omega, ...
                                        k12, k112, k1112, alpha, kappa, ...
                                        theta, gamma_val, h_fun);
    end

    % Same trajectory in z (vectorised x_to_z), for comparison with
    % simulate_case.
    [~, ~, ~, beta, M0] = offhooked_coeffs(d0, d1, l1, l2);
    Y          = (M0 \ x_out.').';
    z_out      = Y;
    z_out(:,5) = Y(:,5) - beta*x_out(:,3).^3;
end

function dx = rhs_x(t, x, zf, d0, d1, l1, l2, omega, k12, k112, k1112, ...
                    alpha, kappa, theta, gamma_val, h_fun)
%RHS_X  Exact kinematics (fields_exact), driven by the controller.
    [f1, f2] = fields_exact(x, d0, d1, l1, l2);
    [u1, u2] = u_ctrl(t, zf, omega, k12, k112, k1112, ...
                      alpha, kappa, theta, gamma_val, h_fun);
    dx = f1*u1 + f2*u2;
end

% ---- CONTROLLER -- u1, u2, amplitudes, worst-case bound ----

function dz = rhs_approx(t, z, zf, g1_fun, g2_fun, omega, k12, k112, k1112, ...
                         alpha, kappa, theta, gamma_val, h_fun)
%RHS_APPROX  Closed-loop right-hand side of the nilpotent approximation,
%   z' = u1*g1(z) + u2*g2(z), with u1, u2 from u_ctrl at the frozen state zf.
%   Used by simulate_case and simulate_case_adaptive.
%
%   g1_fun, g2_fun may be [] -- then the fields come from fields_approx
%   (the normal case). Passing your own handles still works, for
%   experiments with modified fields.

    [u1, u2] = u_ctrl(t, zf, omega, k12, k112, k1112, ...
                      alpha, kappa, theta, gamma_val, h_fun);
    if isempty(g1_fun)
        [g1, g2] = fields_approx(z, alpha, kappa, theta);
    else
        g1 = g1_fun(z, alpha, kappa, theta);
        g2 = g2_fun(z, alpha, kappa, theta);
    end
    dz = g1*u1 + g2*u2;
end

function [u1_val, u2_val] = u_ctrl(t, zf, omega, k12, k112, k1112, ...
                                   alpha, kappa, theta, gamma_val, h_fun)
%U_CTRL  Both control components at time t, for the state zf frozen at the
%   start of the current sampling interval.
%
%   u1 is the driving velocity [m/s], u2 the steering velocity [rad/s].
%   They are PHYSICAL inputs: the same u1, u2 drive the approximate
%   z-system and the exact x-system. Only the STATE needs converting.
%
%   The amplitudes and term sizes come from osc_terms, the same file
%   control_bound uses, so the bound always matches the signal.

    [a1, a2, s12, s112, s1112, T12, T112, T1112] = osc_terms( ...
        zf, alpha, kappa, theta, gamma_val, omega, k12, k112, k1112, h_fun);

    % optional phase offsets (test only): the argument k*omega*t of each
    % oscillating group becomes k*omega*t + ph(i), the same in u1 and u2
    ph = getappdata(0, 'ctrl_phase');
    if isempty(ph), ph = [0 0 0]; end
    th12   = k12*omega*t    + ph(1);
    th112  = k112*omega*t   + ph(2);
    th1112 = k1112*omega*t  + ph(3);

    u1_val = a1 ...
        + T12   * cos(th12) ...
        + T112  * cos(th112) ...
        + 2*T1112 * cos(th1112);

    u2_val = a2 ...
        + T12   * s12   * sin(th12) ...
        - 2*T112  * s112  * cos(2*th112) ...
        - 6*T1112 * s1112 * sin(3*th1112);
end

function [a1, a2, s12, s112, s1112, T12, T112, T1112] = osc_terms( ...
        zf, alpha, kappa, theta, gamma_val, omega, k12, k112, k1112, h_fun)
%OSC_TERMS  Everything the controller needs at the frozen state zf:
%   the constant parts a1, a2, the signs s.. of the bracket amplitudes and
%   the sizes T.. of the three oscillating terms.
%
%   Shared by u_ctrl (the signal) and control_bound (its worst case),
%   so the two can never disagree.

    [a1, a2, a12, a112, a1112] = amps(zf, alpha, kappa, theta, gamma_val, h_fun);

    T12   = (2*k12*omega*abs(a12))^(1/2);
    T112  = (2*k112*omega)^(2/3) * abs(a112)^(1/3);
    T1112 = (k1112*omega)^(3/4) * abs(a1112)^(1/4);

    s12 = sign(a12);  s112 = sign(a112);  s1112 = sign(a1112);
end

function [a1, a2, a12, a112, a1112] = amps(zf, alpha, kappa, theta, gamma_val, h_fun)
    h     = h_fun(zf);
    a1    = -gamma_val*h*zf(1);
    a2    = -gamma_val*h*zf(2);
    a12   = -gamma_val*h*(zf(1)*zf(2) + zf(3));
    a112  = -gamma_val*h*(zf(1)*zf(3) + zf(4) + alpha*zf(2)^3);
    a1112 = -gamma_val*h*( (alpha + 2*kappa - theta)*zf(1)*zf(2)^3 ...
                         + (alpha + 3*kappa)*zf(2)^2*zf(3) ...
                         + zf(1)*zf(4) + zf(5) );
end

function [b1, b2, terms1, terms2] = control_bound(z, alpha, kappa, theta, ...
                                                  gamma_val, omega, k12, k112, k1112, ...
                                                  h_fun)
%CONTROL_BOUND  Worst-case |u1|, |u2| at state z: the constant part plus the
%   magnitude of every oscillating term.
%
%   h_fun is the state-dependent factor of the Lyapunov function (see
%   params.m). It scales ALL five amplitudes, so leaving it out makes this
%   bound wrong by exactly that factor. Omit it only for the plain
%   quadratic V, where h = 1.

    if nargin < 10 || isempty(h_fun)
        h_fun = @(zz) 1;
    end

    [a1, a2, ~, ~, ~, T12, T112, T1112] = osc_terms( ...
        z, alpha, kappa, theta, gamma_val, omega, k12, k112, k1112, h_fun);

    terms1 = [abs(a1), T12,   T112, 2*T1112];
    terms2 = [abs(a2), T12, 2*T112, 6*T1112];

    b1 = sum(terms1);
    b2 = sum(terms2);
end

% ---- GAINS -- admissible gamma and time to converge ----

function g = gamma_max(z, alpha, kappa, theta, omega, k12, k112, k1112, ...
                       u1max, u2max, h_fun)
%GAMMA_MAX  Largest gain for which BOTH velocity limits still hold at state z.
%   Bisection in log space: 200 halvings between 1e-16 and 10.

    if nargin < 11 || isempty(h_fun)
        h_fun = @(zz) 1;
    end

    lo = 1e-16;
    hi = 10;
    for i = 1:200
        mid = sqrt(lo*hi);
        [b1, b2] = control_bound(z, alpha, kappa, theta, mid, omega, ...
                                 k12, k112, k1112, h_fun);
        if b1 <= u1max && b2 <= u2max
            lo = mid;
        else
            hi = mid;
        end
    end
    g = lo;
end

function t = eta_ray(z0, gamma_val, h_fun, fracs)
%ETA_RAY  Time [s] until ||z|| has dropped to frac*||z0||, for any radial V.
%
%   t = eta_ray(P.z0, P.gamma, P.h_fun, [0.1 0.01 0.001])
%
%   For a radial V the closed loop follows (on average over the
%   oscillations) the gradient flow   dr/dt = -gamma*h(r)*r,   r = ||z||,
%   along the straight ray from z0 to the origin. So
%       t(frac) = integral from frac*r0 to r0 of  dr / (gamma*h(r)*r).
%   For the quadratic V (h = 1) this is ln(1/frac)/gamma; for the log and
%   sqrt V it is NOT, which is why a linear estimate 1-eps*gamma per
%   interval is wrong for them.
%
%   (Checked against the full oscillating simulation over 80 s: the two
%   agree to 0.2% of the progress made -- see compare_V.)

    r0   = norm(z0);
    zhat = z0(:) / r0;
    hr   = @(r) arrayfun(@(rr) h_fun(rr*zhat), r);

    t = zeros(size(fracs));
    for i = 1:numel(fracs)
        t(i) = integral(@(r) 1 ./ (gamma_val*hr(r).*r), fracs(i)*r0, r0);
    end
end

% ---- VECTOR FIELDS AND COORDINATES (x_to_z and offhooked_coeffs stay separate files: trailers2.m needs them) ----

function [g1, g2] = fields_approx(z, alpha, kappa, theta)
%FIELDS_APPROX  Vector fields of the nilpotent approximation,
%   z' = u1*g1(z) + u2*g2(z).
%
%   The ONE place these formulas live (calculations, simulate_case,
%   simulate_case_adaptive and the compare_* scripts all use this file).
%   derive_symbolic writes the same fields symbolically.

    g1 = [1; 0; -z(2); -z(3); kappa*z(2)^3 - z(4)];
    g2 = [0; 1; 0; -alpha*z(2)^2; -theta*z(2)*z(3)];
end

function [f1, f2] = fields_exact(x, d0, d1, l1, l2)
%FIELDS_EXACT  Exact kinematics of the car with two off-hooked trailers
%   (paper, eq. (2)-(3)):  x' = u1*f1(x) + u2*f2(x),
%   x = [x1; x2; phi; theta1; theta2], angles in radians.

    x3 = x(3);  x4 = x(4);  x5 = x(5);

    f15 = sin(x3 - x5)/l2 ...
          + (d1 + l1)*sin(x4 - x3)*cos(x4 - x5)/(l1*l2);
    f25 = d0*(d1 + l1)*cos(x3 - x4)*cos(x4 - x5)/(l1*l2) ...
          - d0*cos(x3 - x5)/l2;

    f1 = [cos(x3); sin(x3); 0; sin(x3 - x4)/l1; f15];
    f2 = [0; 0; 1; -d0*cos(x3 - x4)/l1; f25];
end

function x = z_to_x(z, d0, d1, l1, l2)
%Z_TO_X  Privileged coordinates  ->  physical state. Exact inverse of
%   x_to_z.m, useful for plotting where the vehicle actually drives
%   (metres and radians) instead of the abstract z coordinates.
%
%   The inverse is explicit because row 3 of M0 is [0 1 0 0 0], so the
%   heading is simply  phi = x(3) = z(2)  and the cubic correction can be
%   undone without solving anything.

    [~, ~, ~, beta, M0] = offhooked_coeffs(d0, d1, l1, l2);

    z    = z(:);
    y    = z;
    y(5) = z(5) + beta*z(2)^3;      % phi = z(2), see above

    x = M0*y;
end
