function varargout = trailer1(task, varargin)
%TRAILER1  The whole ONE-TRAILER study (car + one off-hooked trailer, r = 3,
%   4 states) in a single file. The numbers are in params.m -- this file
%   has no numbers of its own to edit.
%
%   trailer1                         list of tasks, then run the default one (plots;
%                                    set DEFAULT_TASK below -- also what the Run button does)
%   trailer1('selftest')             does the code still compute what it computed before? (~20 s)
%   trailer1('calculations')         one closed-loop simulation, approximate model, fixed gain
%   trailer1('plots')                states, controls and convergence: approximation vs exact kinematics
%   trailer1('compare_physical')     the same two models as numbers (slow: the full t_end, twice)
%   trailer1('compare_gain')         fixed vs adaptive gain; admissible gain and time vs epsilon
%   trailer1('compare_V')            candidate Lyapunov functions: admissible gain, time, use of the limits
%   trailer1('convergence_map')      which (epsilon, gamma) converge
%   trailer1('search_k', 15)         rank the (k12, k112) pairs
%
%   Every function below can be called the same way from the Command Window,
%   e.g.   [b1, b2] = trailer1('control_bound', P.z0, P.alpha, P.gamma, ...
%                              P.omega, P.k12, P.k112, P.h_fun);
%   (the file is laid out as: tasks, simulators, controller, gains,
%   coordinates and vector fields.)
%
%   The symbolic derivation (Symbolic Math Toolbox) is in derive.m.
%   Do not put this folder and "2 trailers" on the MATLAB path together:
%   both have a params.m.

    DEFAULT_TASK = 'plots';    % what a bare  trailer1  (or the Run button) does;
                               % '' = only list the tasks

    if nargin < 1 || isempty(task)
        show_menu();
        if isempty(DEFAULT_TASK), return, end
        fprintf('Running the default task ''%s'' (change DEFAULT_TASK at the top of trailer1.m)\n\n', DEFAULT_TASK);
        task = DEFAULT_TASK;
    end
    fns = function_table();
    if ~isfield(fns, task)
        error('trailer1:unknownTask', ...
              'Unknown task "%s" -- run trailer1 without arguments to see the list.', task);
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
%   trailer1('name', ...) can reach them (handles to local functions can
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
        'simulate_case', @simulate_case, ...
        'simulate_case_adaptive', @simulate_case_adaptive, ...
        'simulate_case_physical', @simulate_case_physical, ...
        'rhs_approx', @rhs_approx, ...
        'u1_ctrl', @u1_ctrl, ...
        'u2_ctrl', @u2_ctrl, ...
        'osc_terms', @osc_terms, ...
        'amps', @amps, ...
        'control_bound', @control_bound, ...
        'gamma_max', @gamma_max, ...
        'gamma_ray', @gamma_ray, ...
        'eta_ray', @eta_ray, ...
        'fields_approx', @fields_approx, ...
        'fields_exact', @fields_exact, ...
        'x_to_z', @x_to_z, ...
        'z_to_x', @z_to_x, ...
        'heuristic_no_multiple', @heuristic_no_multiple);
end

function show_menu()
    fprintf('\ntrailer1(task, ...)   -- one trailer, r = 3, numbers in params.m\n\n');
    fprintf('  trailer1(''selftest'')           does the code still compute what it did before (~20 s)\n');
    fprintf('  trailer1(''calculations'')       one closed-loop simulation (approximate model)\n');
    fprintf('  trailer1(''plots'')              states / controls / convergence, approx vs exact\n');
    fprintf('  trailer1(''compare_physical'')   approx vs exact kinematics as numbers (slow)\n');
    fprintf('  trailer1(''compare_gain'')       fixed vs adaptive gain, gamma and time vs epsilon\n');
    fprintf('  trailer1(''compare_V'')          candidate Lyapunov functions ranked under the limits\n');
    fprintf('  trailer1(''convergence_map'')    which (epsilon, gamma) converge\n');
    fprintf('  trailer1(''search_k'', 15)       rank the (k12, k112) pairs\n\n');
    fprintf('  Any function in this file can be called as trailer1(''name'', args...).\n');
    fprintf('  Symbolic derivation: derive.m (needs the Symbolic Math Toolbox).\n\n');
end

% ---- TASKS -- what you run ----

function selftest()
%SELFTEST  ONE-TRAILER folder: does it still compute what it computed before?
%
%   Run after ANY edit to the code (about 20 s). Checks
%     1. alpha and z0, and the x -> z -> x round trip,
%     2. the control bound, gamma_max and one value of u1/u2 at z0,
%     3. a 12 s closed-loop run of all three simulators,
%   against reference values, and prints PASSED or which number moved.
%
%   The references were produced on 23 Oct 2026 with params.m as it is
%   now (gamma = 2.1143e-02, eps = 5, k = (3,2), Victoria's log V with
%   c = 10, x0 = (-0.2, -0.2, 45, 45 deg)). Before that they were made for
%   the sqrt V (delta = 0.02, gamma = 3.2028e-03). If you change those ON PURPOSE, this
%   test fails on purpose: check the new numbers, then update the table.

P = params();
n_fail = 0;
K = {P.omega, P.k12, P.k112};

% ---- reference values ----------------------------------------------------
ref.alpha  = 1.132750000000000e-02;
ref.z0     = [-2.000000000000000e-01  7.853981633974483e-01 ...
               2.000000000000000e-01 -2.168101725068444e-02];
ref.bound  = [ 1.853919105523241e-01  3.425117712486228e-01];
ref.gmax   =   3.388979397850971e-02;
ref.uctrl  = [ 9.867208209915458e-02  5.775493891793783e-02];   % at t = 0.37 s

% 12 s runs; the ODE solver may differ in the last digits between MATLAB
% and Octave, hence the looser tolerance for these
ref.z_approx = [-2.447674483114e-01  7.909991514118e-01  2.357674500212e-01 -1.017954713026e-02];
ref.u_approx = [ 1.874920948451e-01  3.437566059164e-01];
ref.z_adapt  = [-2.478831391779e-01  7.804096868508e-01  2.376370100693e-01 -8.390606257002e-03];
ref.u_adapt  = [ 2.200000000000e-01  4.045503847606e-01];
ref.x_exact  = [-2.310717573053e-01 -2.319120710183e-01  7.909263190036e-01 -3.119159139290e-01];
ref.u_exact  = [ 1.864540281977e-01  3.428780580405e-01];

% ---- 1. geometry-only numbers and the coordinate map ----------------------
n_fail = n_fail + selftest_check('alpha', P.alpha, ref.alpha, 1e-12);
n_fail = n_fail + selftest_check('z0 = x0_to_z0(x0)', P.z0(:).', ref.z0, 1e-12);

xr = [P.x0(1); P.x0(2); deg2rad(P.x0(3)); deg2rad(P.x0(4))];
n_fail = n_fail + selftest_check('x -> z -> x round trip', ...
    z_to_x(x_to_z(xr, P.d0, P.l1), P.d0, P.l1).', xr.', 1e-12);
n_fail = n_fail + selftest_check('x_to_z = x0_to_z0', x_to_z(xr, P.d0, P.l1).', P.z0(:).', 1e-12);

% ---- 2. controller at z0 ---------------------------------------------------
[b1, b2] = control_bound(P.z0, P.alpha, P.gamma, K{:}, P.h_fun);
n_fail = n_fail + selftest_check('control_bound at z0', [b1 b2], ref.bound, 1e-12);
g = gamma_max(P.z0, P.alpha, K{:}, P.u1_max, P.u2_max, P.h_fun);
n_fail = n_fail + selftest_check('gamma_max at z0', g, ref.gmax, 1e-10);
u1 = u1_ctrl(0.37, P.z0, K{:}, P.alpha, P.gamma, P.h_fun);
u2 = u2_ctrl(0.37, P.z0, K{:}, P.alpha, P.gamma, P.h_fun);
n_fail = n_fail + selftest_check('u1/u2 at t = 0.37 s', [u1 u2], ref.uctrl, 1e-12);

% ---- 3. 12 s closed loop, all three simulators -------------------------------
T = {0, 12, P.dt};                         % t_start, t_end, dt
fprintf('\nsimulating 12 s (approx / adaptive / exact) ...\n');

[~, z, u1, u2] = simulate_case(P.z0, P.alpha, P.gamma, P.epsilon, K{:}, T{:}, P.h_fun);
n_fail = n_fail + selftest_check('approx   z(12 s)', z(end,:), ref.z_approx, 1e-7);
n_fail = n_fail + selftest_check('approx   max|u|', [max(abs(u1)) max(abs(u2))], ref.u_approx, 1e-7);

[~, z, u1, u2] = simulate_case_adaptive(P.z0, P.alpha, P.epsilon, K{:}, T{:}, P.h_fun, ...
                                        [], [], P.u1_max, P.u2_max, P.cap_product);
n_fail = n_fail + selftest_check('adaptive z(12 s)', z(end,:), ref.z_adapt, 1e-7);
n_fail = n_fail + selftest_check('adaptive max|u|', [max(abs(u1)) max(abs(u2))], ref.u_adapt, 1e-7);

[~, x, ~, u1, u2] = simulate_case_physical(xr, P.d0, P.l1, P.alpha, P.gamma, ...
                                           P.epsilon, K{:}, T{:}, P.h_fun);
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
% CALCULATIONS  One closed-loop simulation of the ONE-TRAILER system
% (r = 3, 4 states) with the numbers from params.m.
%
%   Runs the approximate (nilpotent) model at the fixed gain. As a command
%   it leaves t_sim, z_sim, u1_sim, u2_sim, normz_sim in the base workspace;
%   with outputs, [t,z,u1,u2,nz] = trailer1('calculations'), it returns them.
%   Uses the Lyapunov function set in params.m (the gain there is tuned for
%   it). No Symbolic Toolbox needed: the symbolic derivation is in derive.m.

    P = params();
    d0 = P.d0;  l1 = P.l1;
    alpha = P.alpha;
    gamma = P.gamma;  epsilon = P.epsilon;  omega = P.omega;
    k12 = P.k12;  k112 = P.k112;
    z0 = P.z0;  t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;
    U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

    fprintf('geometry : d0=%.3f  l1=%.3f  [m]\n', d0, l1);
    fprintf('coeff    : alpha=%.6f\n', alpha);

    % ---- One closed-loop simulation (approximate model, fixed gain, V from params.m) ----
    h_fun = P.h_fun;

    fprintf('\n========================================\n')
    fprintf('  Simulation\n')
    fprintf('  z0      = [%.4f; %.4f; %.4f; %.4f],  ||z0|| = %.4f\n', z0, norm(z0));
    fprintf('  time    : [%.1f, %.1f], dt = %.1e\n', t_start, t_end, dt);
    fprintf('  epsilon = %.2f, omega = %.4f, gamma = %.4e\n', epsilon, omega, gamma);
    fprintf('  V       : %s\n', P.lyap_name);
    fprintf('  k       = (%d, %d)   [k12, k112 -- no k1112 for r=3]\n', k12, k112);
    fprintf('========================================\n')

    [t_sim, z_sim, u1_sim, u2_sim] = simulate_case( ...
        z0, alpha, gamma, ...
        epsilon, omega, k12, k112, ...
        t_start, t_end, dt, h_fun);

    fprintf('  %d steps\n', length(t_sim));
    fprintf('  final state : [%.6f; %.6f; %.6f; %.6f]\n', z_sim(end,:));
    fprintf('  final norm  : %.6f\n', norm(z_sim(end,:)));
    fprintf('  max |u1|    : %.3f  (limit %.2f)\n', max(abs(u1_sim)), U1_MAX);
    fprintf('  max |u2|    : %.3f  (limit %.2f)\n', max(abs(u2_sim)), U2_MAX);

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
% PLOTS  Figures for the ONE-TRAILER (r=3, 4-state) simulation.
%
%   Draws the approximate model (calculations, run for you) and
%   the SAME controller on the exact kinematics, with the same Lyapunov
%   function from params.m in both runs: states, controls against the
%   TurtleBot3 limits, convergence of ||z||.
%
%   Run:  trailer1('plots')

    P = params();
    U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

    col     = [0.00 0.20 0.70];
    col_ex  = [0.85 0.33 0.10];

    [t_sim, z_sim, u1_sim, u2_sim, normz_sim] = calculations();

    % The SAME controller on the exact kinematics. u1, u2 are not "converted"
    % here: the formula is identical, it is just fed the state of a different
    % plant, so the commanded values come out different. This is what the real
    % robot would actually be asked to do.
    fprintf('(running the exact-kinematics model for comparison ...)\n');
    x0_rad = [P.x0(1); P.x0(2); deg2rad(P.x0(3)); deg2rad(P.x0(4))];
    [t_ex, ~, z_ex, u1_ex, u2_ex] = simulate_case_physical( ...
        x0_rad, P.d0, P.l1, P.alpha, P.gamma, P.epsilon, P.omega, ...
        P.k12, P.k112, P.t_start, P.t_end, P.dt, P.h_fun);

    % ---------------- state trajectories ----------------------------------------
    figure('Name', 'State trajectories (1 trailer)', 'NumberTitle', 'off');
    labels = {'z_1', 'z_2', 'z_3', 'z_4'};
    for i = 1:4
        subplot(4,1,i)
        plot(t_sim, z_sim(:,i), 'Color', col, 'LineWidth', 1.2); hold on
        plot(t_ex,  z_ex(:,i),  '--', 'Color', col_ex, 'LineWidth', 1.0);
        grid on; ylabel(labels{i});
        if i == 1
            title('State trajectories (solid: approximation, dashed: exact kinematics)');
        end
        if i == 4
            xlabel('t [s]');
        end
    end

    % ---------------- control inputs ---------------------------------------------
    figure('Name', 'Control inputs (1 trailer)', 'NumberTitle', 'off');
    subplot(2,1,1)
    plot(t_sim, u1_sim, 'Color', col,    'LineWidth', 1.0); hold on
    plot(t_ex,  u1_ex,  'Color', col_ex, 'LineWidth', 1.0);
    yline(U1_MAX, 'r--'); yline(-U1_MAX, 'r--');
    grid on; ylabel('u_1 [m/s]');
    title('Control inputs (blue: approximation, orange: exact; red dashed: limits)');

    subplot(2,1,2)
    plot(t_sim, u2_sim, 'Color', col,    'LineWidth', 1.0); hold on
    plot(t_ex,  u2_ex,  'Color', col_ex, 'LineWidth', 1.0);
    yline(U2_MAX, 'r--'); yline(-U2_MAX, 'r--');
    grid on; xlabel('t [s]'); ylabel('u_2 [rad/s]');

    % ---------------- convergence -------------------------------------------------
    figure('Name', 'Convergence (1 trailer)', 'NumberTitle', 'off');
    semilogy(t_sim, max(normz_sim, realmin), 'Color', col, 'LineWidth', 1.5); hold on
    semilogy(t_ex, max(vecnorm(z_ex, 2, 2), realmin), 'Color', col_ex, 'LineWidth', 1.5);
    grid on; xlabel('t [s]'); ylabel('||z(t)||');
    title('Convergence');
    legend('approximation', 'exact kinematics', 'Location', 'best');

    r1 = max(abs(u1_sim))/U1_MAX;
    r2 = max(abs(u2_sim))/U2_MAX;
    r1e = max(abs(u1_ex))/U1_MAX;
    r2e = max(abs(u2_ex))/U2_MAX;
    fprintf('\n  V : %s\n', P.lyap_name);
    fprintf('\n  max |u1| = %8.4f m/s    limit %.2f    -> %6.1f %% of the limit\n', max(abs(u1_sim)), U1_MAX, 100*r1);
    fprintf('  max |u2| = %8.4f rad/s  limit %.2f    -> %6.1f %% of the limit\n', max(abs(u2_sim)), U2_MAX, 100*r2);
    fprintf('\n  --- same controller, EXACT kinematics (what the robot does) ---\n');
    fprintf('  max |u1| = %8.4f m/s    limit %.2f    -> %6.1f %% of the limit\n', max(abs(u1_ex)), U1_MAX, 100*r1e);
    fprintf('  max |u2| = %8.4f rad/s  limit %.2f    -> %6.1f %% of the limit\n', max(abs(u2_ex)), U2_MAX, 100*r2e);
    fprintf('  final ||z|| : approximation %.3e,  exact %.3e\n', ...
            norm(z_sim(end,:)), norm(z_ex(end,:)));
    if r1e > 1 + 1e-4 || r2e > 1 + 1e-4      % 1e-4: P.gamma is rounded to 5 digits
        fprintf(2, '  *** the EXACT model exceeds a limit -- see OBALKA_1PRIVES.md ***\n');
    end
end

function compare_physical()
% COMPARE_PHYSICAL
% Same controller, two plants, ONE TRAILER (r = 3):
%   (a) the nilpotent approximation  z' = u1*g1(z) + u2*g2(z)
%       -- what simulate_case / calculations show
%   (b) the exact kinematics         x' = u1*f1(x) + u2*f2(x)
%       -- what the real robot, and Gazebo, are
% The gap is the remainder the approximation drops. Run this BEFORE
% trusting any retune and before porting anything to Gazebo: at the
% two-trailer case a setting that looked perfect on (a) diverged on (b).

    P = params();
    d0 = P.d0;  l1 = P.l1;  alpha = P.alpha;
    gamma = P.gamma;  epsilon = P.epsilon;  omega = P.omega;
    k12 = P.k12;  k112 = P.k112;
    z0 = P.z0;  t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;

    h_fun  = P.h_fun;

    % CAREFUL: P.x0 stores the angles in DEGREES (params.m converts them), but
    % the physical simulation integrates radians.
    x0 = [P.x0(1); P.x0(2); deg2rad(P.x0(3)); deg2rad(P.x0(4))];

    fprintf('========================================\n');
    fprintf('  x0      = [%.3f; %.3f; %.1f deg; %.1f deg]\n', ...
            P.x0(1), P.x0(2), P.x0(3), P.x0(4));
    fprintf('  z0      = [%.4f; %.4f; %.4f; %.4f],  ||z0|| = %.4f\n', z0, norm(z0));
    fprintf('  gamma   = %.4e,  eps = %.2f,  k = (%d,%d)\n', gamma, epsilon, k12, k112);
    fprintf('  V       : %s\n', P.lyap_name);
    fprintf('  time    : [%.0f, %.0f]\n', t_start, t_end);
    fprintf('========================================\n');

    fprintf('\nrunning approximate (z) model ...\n');
    [t_a, z_a, u1_a, u2_a] = simulate_case( ...
        z0, alpha, gamma, epsilon, omega, k12, k112, ...
        t_start, t_end, dt, h_fun);

    fprintf('running exact (x) model ...\n');
    [t_p, x_p, z_p, u1_p, u2_p] = simulate_case_physical( ...
        x0, d0, l1, alpha, gamma, epsilon, omega, k12, k112, ...
        t_start, t_end, dt, h_fun);

    % ---------------------------------------------------------------- results
    na = sqrt(sum(z_a.^2, 2));   np_ = sqrt(sum(z_p.^2, 2));
    hit = @(t, n, f) deal_hit(t, n, f*norm(z0));

    fprintf('\n--- convergence ---\n');
    for f = [0.1 0.01 0.001]
        fprintf('  ||z|| -> %6.1f%% of ||z0||   approx %10s   exact %10s\n', ...
                f*100, hit(t_a, na, f), hit(t_p, np_, f));
    end
    fprintf('  final ||z||                  approx %10.3e   exact %10.3e\n', ...
            na(end), np_(end));

    fprintf('\n--- final physical state (exact model) ---\n');
    fprintf('  pos (%+.4f, %+.4f) m   phi %+.2f deg   theta1 %+.2f deg\n', ...
            x_p(end,1), x_p(end,2), rad2deg(x_p(end,3)), rad2deg(x_p(end,4)));

    fprintf('\n--- controls (same formula in both runs) ---\n');
    fprintf('  approx : max|u1| = %.4f,  max|u2| = %.4f\n', max(abs(u1_a)), max(abs(u2_a)));
    fprintf('  exact  : max|u1| = %.4f,  max|u2| = %.4f\n', max(abs(u1_p)), max(abs(u2_p)));
    fprintf('  limits : u1_max  = %.4f,  u2_max  = %.4f\n', P.u1_max, P.u2_max);
    if max(abs(u1_p)) > P.u1_max || max(abs(u2_p)) > P.u2_max
        fprintf(2, '  *** the EXACT model exceeds a limit -- do not trust this tuning ***\n');
    end

    t_c  = linspace(t_start, min(t_a(end), t_p(end)), 2000).';
    drift = sqrt(sum((interp1(t_a, z_a, t_c) - interp1(t_p, z_p, t_c)).^2, 2));
    fprintf('\n--- cost of the nilpotent approximation ---\n');
    fprintf('  max ||z_exact - z_approx|| = %.5f  (%.1f%% of ||z0||), at t = %.1f s\n', ...
            max(drift), 100*max(drift)/norm(z0), t_c(find(drift == max(drift), 1)));

    % ------------------------------------------------------------------ plots
    figure('Name', 'one trailer: approximate vs exact', 'Color', 'w');

    subplot(2,2,1);
    semilogy(t_a, na, 'LineWidth', 1.1); hold on
    semilogy(t_p, np_, 'LineWidth', 1.1);
    grid on; xlabel('t [s]'); ylabel('||z||');
    legend('approximate (z model)', 'exact (x model)', 'Location', 'best');
    title('norm in privileged coordinates');

    subplot(2,2,2);
    plot(x_p(:,1), x_p(:,2), 'LineWidth', 1.1); hold on
    plot(x_p(1,1), x_p(1,2), 'o'); plot(0, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);
    grid on; axis equal; xlabel('x_1 [m]'); ylabel('x_2 [m]');
    legend('path', 'start', 'goal', 'Location', 'best');
    title('where the vehicle actually drives');

    subplot(2,2,3);
    plot(t_p, rad2deg(x_p(:,3:4)), 'LineWidth', 1.1);
    grid on; xlabel('t [s]'); ylabel('angle [deg]');
    legend('\phi', '\theta_1', 'Location', 'best');
    title('physical angles, exact model');

    subplot(2,2,4);
    plot(t_p, u1_p, 'LineWidth', 0.8); hold on
    yline( P.u1_max, '--'); yline(-P.u1_max, '--');
    grid on; xlabel('t [s]'); ylabel('u_1 [m/s]');
    title('driving velocity vs limit');

    fprintf('\n=== Done. ===\n');
end

function s = deal_hit(t, n, target)
    i = find(n <= target, 1);
    if isempty(i), s = 'never'; else, s = sprintf('%.0f s', t(i)); end
end

function compare_gain()
% COMPARE_GAIN  Fixed vs. adaptive, hardware-aware gain, and how the
% admissible gain depends on epsilon. ONE-TRAILER (r=3, 4-state) version --
% see "2 trailers/compare_gain.m" for the two-trailer version.
%
%   fixed    : gamma = P.gamma for the whole run.
%   adaptive : gamma recomputed every interval = largest value that still
%              respects the robot's |u1|,|u2| limits (gamma_max),
%              additionally capped so eps*gamma*h(z) <= P.cap_product.
%
% Uses the Lyapunov function from params.m. Time-to-converge estimates use
% eta_ray (gradient flow along the ray, valid for ANY radial V); the old
% estimate (1-eps*gamma)^N assumed the quadratic V.
%
% Run:  trailer1('compare_gain')

    clc; close all;

    P = params();
    alpha = P.alpha;
    k12 = P.k12;  k112 = P.k112;
    z0 = P.z0;  t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;
    U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

    h_fun  = P.h_fun;

    fprintf('geometry : d0=%.3f  l1=%.3f  [m]\n', P.d0, P.l1);
    fprintf('z0 = [%g %g %g %g],  eps = %g,  k = (%d,%d)\n', z0, P.epsilon, k12, k112);
    fprintf('V  : %s\n', P.lyap_name);

    % ---------------- fixed gain ------------------------------------------------
    fprintf('\n--- fixed gain, gamma = %g ---\n', P.gamma);
    [t_f, z_f, u1_f, u2_f] = simulate_case(z0, alpha, P.gamma, ...
        P.epsilon, P.omega, k12, k112, t_start, t_end, dt, h_fun);
    fprintf('  final norm : %.6f\n', norm(z_f(end,:)));
    fprintf('  max |u1|   : %.3f  (limit %.2f)\n', max(abs(u1_f)), U1_MAX);
    fprintf('  max |u2|   : %.3f  (limit %.2f)\n', max(abs(u2_f)), U2_MAX);

    % ---------------- adaptive gain ---------------------------------------------
    fprintf('\n--- adaptive gain, cap eps*gamma*h <= %g ---\n', P.cap_product);
    [t_a, z_a, u1_a, u2_a, gamma_a] = simulate_case_adaptive(z0, alpha, ...
        P.epsilon, P.omega, k12, k112, t_start, t_end, dt, h_fun, [], [], ...
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
        fprintf('  ||z|| down to %5.1f%% of ||z0|| : %8.0f s = %6.1f min\n', ...
                targets(i)*100, t_need(i), t_need(i)/60);
    end

    col_f = [0.00 0.20 0.70];
    col_a = [0.85 0.33 0.10];

    % ---------------- gamma_max and time-to-converge vs epsilon -----------------
    eps_sweep    = [0.5 1 2 3 4 5 6 8 12 16 32];
    n_sweep      = numel(eps_sweep);
    gmax_sweep   = zeros(1, n_sweep);
    min10_sweep  = nan(1, n_sweep);
    min1_sweep   = nan(1, n_sweep);

    fprintf('\n===== gamma_max and time-to-converge vs epsilon =====\n');
    fprintf('(gamma = largest gain admissible along the whole ray z0 -> 0 for the velocity\n');
    fprintf(' limits AND eps*gamma*h <= %.2f, see gamma_ray; time from eta_ray)\n', P.cap_product);
    fprintf('%10s %14s %14s %14s %14s\n', 'eps [s]', 'gamma', 'eps*g*hmax', 'min to 10%', 'min to 1%');
    for i = 1:n_sweep
        Q = params(eps_sweep(i));
        [gmax_sweep(i), hmax] = gamma_ray(z0, alpha, Q.omega, Q.epsilon, k12, k112, ...
                                          U1_MAX, U2_MAX, h_fun, P.cap_product);
        t = eta_ray(z0, gmax_sweep(i), h_fun, [0.1 0.01]) / 60;
        min10_sweep(i) = t(1);
        min1_sweep(i)  = t(2);
        fprintf('%10.2f %14.3e %14.3e %14.1f %14.1f\n', Q.epsilon, gmax_sweep(i), ...
                Q.epsilon*gmax_sweep(i)*hmax, min10_sweep(i), min1_sweep(i));
    end
    fprintf(['NOTE: the sweep uses the approximate model. params.m records eps = 2..7\n' ...
             '      as verified against the exact kinematics for this start state; for\n' ...
             '      larger eps (and for two trailers beyond eps ~ 4 s the approximation\n' ...
             '      broke down) check with compare_physical before believing it.\n']);

    figure('Name', 'gamma_max and ETA vs epsilon (1 trailer)', 'NumberTitle', 'off');
    subplot(1,2,1)
    loglog(eps_sweep, gmax_sweep, 'o-', 'Color', col_a, 'LineWidth', 1.5);
    grid on; xlabel('\epsilon  [s]'); ylabel('\gamma  (admissible along the ray)');
    title('Largest admissible gain vs. period');

    subplot(1,2,2)
    loglog(eps_sweep, min10_sweep, 'o-',  'Color', col_a, 'LineWidth', 1.5); hold on
    loglog(eps_sweep, min1_sweep,  's--', 'Color', col_a, 'LineWidth', 1.0);
    grid on; xlabel('\epsilon  [s]'); ylabel('minutes to converge');
    legend('to 10% of ||z_0||', 'to 1% of ||z_0||', 'Location', 'best');
    title('Estimated time to converge at the admissible gain');

    % ---------------- plot: ||z(t)|| --------------------------------------------
    figure('Name', 'Fixed vs adaptive gain -- convergence (1 trailer)', 'NumberTitle', 'off');
    semilogy(t_f, vecnorm(z_f, 2, 2), 'Color', col_f, 'LineWidth', 1.5); hold on
    semilogy(t_a, vecnorm(z_a, 2, 2), 'Color', col_a, 'LineWidth', 1.5);
    grid on; xlabel('t [s]'); ylabel('||z(t)||');
    legend('fixed gain', 'adaptive gain', 'Location', 'best');
    title('Convergence: fixed vs. adaptive gain');

    % ---------------- plot: controls, both cases --------------------------------
    figure('Name', 'Fixed vs adaptive gain -- controls (1 trailer)', 'NumberTitle', 'off');
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
    figure('Name', 'Adaptive gain over time (1 trailer)', 'NumberTitle', 'off');
    plot(t_a, gamma_a, 'Color', col_a, 'LineWidth', 1.5); hold on
    yline(P.gamma, 'b--');
    grid on; xlabel('t [s]'); ylabel('\gamma used');
    legend('adaptive \gamma(t)', 'fixed \gamma', 'Location', 'best');
    title('Gain chosen every interval (hardware-limited, capped \epsilon\gamma h <= cap)');

    fprintf('\n=== Done. ===\n');
end

function compare_V()
% COMPARE_V  Candidate Lyapunov functions, ranked under the robot's limits.
% ONE-TRAILER (r=3, 4-state) version -- see "2 trailers/compare_V.m".
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
%      can be found by walking that ray (gamma_ray).
%   2. The closed loop approximates the gradient flow  dr/dt = -gamma*h(r)*r,
%      a 1-D ODE, so the convergence time follows without simulating the
%      oscillation (eta_ray). These times are therefore ESTIMATES; the row
%      "in use" is what trailer1('compare_physical') checks by simulation.
%
% Per candidate, gamma is the largest value that keeps
%   (1) the control bound within the limits at EVERY point of the ray, and
%   (2) eps*gamma*h <= P.cap_product on the ray (contraction, see params.m).
% The row marked "in use" is whatever params.m currently sets.
%
% Run:  trailer1('compare_V')

    P = params();
    alpha = P.alpha;  omega = P.omega;  k12 = P.k12;  k112 = P.k112;
    U1 = P.u1_max;  U2 = P.u2_max;
    z0 = P.z0;  r0 = norm(z0);  zhat = z0/r0;

    CAND = { ...
        ['in use: ' P.lyap_name],                 P.h_fun; ...
        'quadratic  V = gamma*r^2/2',             @(z) 1; ...
        'log c=1    V = gamma/(2c)log(1+cr^2)',   @(z) 1/(1 +   1*norm(z)^2); ...
        'log c=10   (Grushkovskaya''s start)',    @(z) 1/(1 +  10*norm(z)^2); ...
        'log c=100',                              @(z) 1/(1 + 100*norm(z)^2); ...
        'sqrt d=0.05 V = gamma(sqrt(r^2+d^2)-d)', @(z) 1/sqrt(norm(z)^2 + 0.05^2); ...
        'sqrt d=0.02',                            @(z) 1/sqrt(norm(z)^2 + 0.02^2); ...
        'sqrt d=0.01',                            @(z) 1/sqrt(norm(z)^2 + 0.01^2); ...
        'sqrt d=0.005',                           @(z) 1/sqrt(norm(z)^2 + 0.005^2) };

    rs = r0 * logspace(-3, 0, 200);

    fprintf('||z0|| = %.5f,  limits |u1|<=%.2f, |u2|<=%.2f,  k=(%d,%d), eps=%.1f, cap=%.2f\n\n', ...
            r0, U1, U2, k12, k112, P.epsilon, P.cap_product);
    nC  = size(CAND, 1);
    res = zeros(nC, 7);     % gamma, t10, t1, t01, eps*g*hmax, max u1 use, max u2 use
    for c = 1:nC
        h = CAND{c, 2};

        % (1)+(2): largest gain admissible at EVERY point of the ray
        [g, hmax] = gamma_ray(z0, alpha, omega, P.epsilon, k12, k112, ...
                              U1, U2, h, P.cap_product);

        % how much of each limit the chosen gain uses (worst point of the ray)
        use1 = 0;  use2 = 0;
        for r = rs
            [b1, b2] = control_bound(r*zhat, alpha, g, omega, k12, k112, h);
            use1 = max(use1, b1/U1);  use2 = max(use2, b2/U2);
        end

        t = eta_ray(z0, g, h, [0.10 0.01 0.001]) / 60;      % minutes
        res(c, :) = [g, t, P.epsilon*g*hmax, use1, use2];
    end
    base = res(2, 4);                % row 2 = quadratic V, time to 0.1%

    fprintf('%-40s %10s %8s %8s %8s %8s %9s %6s %6s\n', 'Lyapunov function', 'gamma', ...
            't->10%', 't->1%', 't->0.1%', 'speedup', 'eps*g*h', 'u1 %', 'u2 %');
    for c = 1:nC
        fprintf('%-40s %10.3e %5.1f m %5.1f m %5.1f m %7.2fx %9.2e %5.0f%% %5.0f%%\n', ...
                CAND{c,1}, res(c,1), res(c,2), res(c,3), res(c,4), base/res(c,4), ...
                res(c,5), 100*res(c,6), 100*res(c,7));
    end
    fprintf(['\nt = estimated minutes until ||z|| has dropped to 10%% / 1%% / 0.1%% of ||z0||.\n' ...
             'speedup = time to 0.1%% of the quadratic V divided by this row''s.\n' ...
             'eps*g*h = largest eps*gamma*h on the ray (contraction, must stay <= cap).\n' ...
             'u1 %%, u2 %% = largest share of each velocity limit used on the ray;\n' ...
             '100%% marks the limit that sets gamma (or the cap does, if eps*g*h = cap).\n']);
end

function convergence_map()
% CONVERGENCE_MAP  Does (epsilon, gamma) converge, overshoot, or diverge?
% ONE-TRAILER (r=3, 4-state) version -- see "2 trailers/convergence_map.m"
% for the full explanation of the method and the linear-bound reference
% lines. Same idea, just the 4-state system and 2 multipliers (k12,k112).
% Uses the Lyapunov function from params.m, so the map is for the
% controller actually in use.
%
% Run:  trailer1('convergence_map')

    clc; close all;

    P = params();
    alpha = P.alpha;
    k12 = P.k12;  k112 = P.k112;
    z0 = P.z0;  t_start = P.t_start;

    h_fun  = P.h_fun;

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
        Q       = params(eps_list(i));
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
end

function T = search_k(kmax, eps_override)
%SEARCH_K  Heuristically-admissible (k12,k112) pairs, ranked by admissible
%   gain. ONE-TRAILER (r=3) version -- mirrors "2 trailers/search_k.m",
%   but with only two multipliers and a HEURISTIC screen instead of the
%   paper's proven condition C1 (see heuristic_no_multiple: this is NOT
%   confirmed to be the right condition for this reduced case).
%
%   T = search_k(15)        search k = 1..15, epsilon from params.m
%   T = search_k(15, 8.0)   with epsilon = 8 s

    if nargin < 1, kmax = 15; end
    P = params();
    if nargin >= 2
        P = params(eps_override);
    end
    alpha = P.alpha;
    gamma = P.gamma;  epsilon = P.epsilon;  w = P.omega;
    z0 = P.z0;
    h_fun  = P.h_fun;
    U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

    fprintf('geometry : d0=%.3f  l1=%.3f  [m]\n', P.d0, P.l1);
    fprintf('coeff    : alpha=%.6f\n', alpha);
    fprintf('V        : %s,  max_u1 column at gamma = %.4e\n', P.lyap_name, gamma);

    rows = zeros(0, 4);
    for a = 1:kmax
        for b = 1:kmax
            if heuristic_no_multiple(a, b)
                g  = gamma_max(z0, alpha, w, a, b, U1_MAX, U2_MAX, h_fun);
                u1 = control_bound(z0, alpha, gamma, w, a, b, h_fun);
                rows(end+1, :) = [a, b, g, u1];   %#ok<AGROW>
            end
        end
    end

    [~, idx] = sort(rows(:, 3), 'descend');
    rows = rows(idx, :);

    T = array2table(rows, 'VariableNames', ...
        {'k12', 'k112', 'gamma_max', 'max_u1'});

    base = gamma_max(z0, alpha, w, P.k12, P.k112, U1_MAX, U2_MAX, h_fun);
    fprintf('\n %d heuristic-admissible pairs with k <= %d,  eps = %.2f s\n', ...
            size(rows, 1), kmax, epsilon);
    fprintf(' current default (%d,%d) gives gamma_max = %.3e\n\n', ...
            P.k12, P.k112, base);
    fprintf('%6s %6s %14s %12s %10s\n', ...
            'k12', 'k112', 'gamma_max', 'max|u1|', 'vs default');
    for i = 1:min(10, size(rows, 1))
        fprintf('%6d %6d %14.4e %12.3f %9.2fx\n', ...
                rows(i,1), rows(i,2), rows(i,3), rows(i,4), rows(i,3)/base);
    end
    fprintf('\n(condition used is a heuristic, not the proven C1 -- see\n');
    fprintf('heuristic_no_multiple and the note in params.m)\n\n');
end

% ---- SIMULATORS -- the three closed loops (approximate / adaptive gain / exact kinematics) ----

function [t_out, z_out, u1_out, u2_out] = simulate_case( ...
        z0, alpha, gamma_val, ...
        epsilon, omega, k12, k112, ...
        t_start, t_end, dt, h_fun, g1_fun, g2_fun)
%SIMULATE_CASE  Sampling scheme: the amplitudes are computed from the state
%   at the start of each interval and held frozen for its whole length.
%   ONE-TRAILER (r=3), nilpotent approximation  z' = u1*g1(z) + u2*g2(z),
%   fixed gain gamma_val.
%
%   g1_fun, g2_fun are optional: leave them out (or pass []) and the fields
%   come from fields_approx.
%   For the exact kinematics use simulate_case_physical, for the
%   state-dependent gain simulate_case_adaptive.

    if nargin < 12, g1_fun = []; end
    if nargin < 13, g2_fun = []; end

    % the -1e-9 stops float noise (e.g. 3.0000000000000004) adding an empty interval
    N_int   = ceil((t_end - t_start)/epsilon - 1e-9);
    t_out   = [];
    z_out   = [];
    z_cur   = z0(:);
    frozen  = zeros(4, N_int);
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', dt);

    for j = 1:N_int
        t_j   = t_start + (j-1)*epsilon;
        t_jp1 = min(t_j + epsilon, t_end);
        if t_j >= t_end
            frozen  = frozen(:, 1:j-1);
            break
        end

        z_frozen     = z_cur;
        frozen(:, j) = z_frozen;

        odefun = @(t, z) rhs_approx(t, z, z_frozen, g1_fun, g2_fun, omega, ...
                                    k12, k112, alpha, gamma_val, h_fun);

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
        u1_out(i) = u1_ctrl(t_out(i), zf, omega, k12, k112, alpha, gamma_val, h_fun);
        u2_out(i) = u2_ctrl(t_out(i), zf, omega, k12, k112, alpha, gamma_val, h_fun);
    end
end

function [t_out, z_out, u1_out, u2_out, gamma_out] = simulate_case_adaptive( ...
        z0, alpha, ...
        epsilon, omega, k12, k112, ...
        t_start, t_end, dt, h_fun, g1_fun, g2_fun, ...
        u1_max, u2_max, cap_product)
%SIMULATE_CASE_ADAPTIVE  ONE-TRAILER (r=3) version. Same idea as the
%   2-trailer version: gamma is recomputed at the start of every interval
%   as the largest value that (a) keeps the commanded |u1|,|u2| within the
%   robot's hardware limits at the frozen state (gamma_max), and (b)
%   keeps eps*gamma*h(z) under cap_product.
%
%   g1_fun, g2_fun may be [] -- then the fields come from fields_approx.
%   gamma_out is the gain in use at every output time.

    % the -1e-9 stops float noise (e.g. 3.0000000000000004) adding an empty interval
    N_int     = ceil((t_end - t_start)/epsilon - 1e-9);
    t_out     = [];
    z_out     = [];
    z_cur     = z0(:);
    frozen    = zeros(4, N_int);
    gamma_int = zeros(1, N_int);
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

        % h_fun must go in here too: the gain is chosen from the same
        % amplitudes the controller then uses, so a mismatch would pick a
        % gain that violates the limits by exactly the factor h.
        g_hw = gamma_max(z_frozen, alpha, omega, k12, k112, u1_max, u2_max, h_fun);
        % contraction cap eps*gamma*h <= cap_product (params.m, rule (2)).
        % For the quadratic V h = 1 and this is the old eps*gamma <= cap;
        % for the sqrt / log V h varies with the state and must be included.
        g_j  = min(g_hw, cap_product/(epsilon*h_fun(z_frozen)));
        gamma_int(j) = g_j;

        odefun = @(t, z) rhs_approx(t, z, z_frozen, g1_fun, g2_fun, omega, ...
                                    k12, k112, alpha, g_j, h_fun);

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
        u1_out(i) = u1_ctrl(t_out(i), zf, omega, k12, k112, alpha, g_j, h_fun);
        u2_out(i) = u2_ctrl(t_out(i), zf, omega, k12, k112, alpha, g_j, h_fun);
    end
end

function [t_out, x_out, z_out, u1_out, u2_out] = simulate_case_physical( ...
        x0, d0, l1, alpha, gamma_val, epsilon, omega, k12, k112, ...
        t_start, t_end, dt, h_fun)
%SIMULATE_CASE_PHYSICAL  Closed loop on the REAL model. ONE-TRAILER (r=3).
%
%   simulate_case integrates the nilpotent approximation
%       z' = u1*g1(z) + u2*g2(z),
%   i.e. the truncated system the controller was designed on. This file
%   integrates the exact kinematics instead,
%       x' = u1*f1(x) + u2*f2(x),
%   driven by exactly the same controller. The gap between the two runs is
%   the remainder the approximation throws away -- which is what the real
%   robot, and Gazebo, will actually feel.
%
%   x0 = [x1; x2; phi; theta1] in RADIANS.
%
%   Loop per sampling interval (also the loop the robot has to run):
%     1. read the physical state x
%     2. z = x_to_z(x)              <- state conversion
%     3. u1, u2 from the controller at that frozen z
%     4. apply u1, u2               <- NO conversion back; u1, u2 are
%                                      already the physical velocities

    % the -1e-9 stops float noise (e.g. 3.0000000000000004) adding an empty interval
    N_int   = ceil((t_end - t_start)/epsilon - 1e-9);
    t_out   = [];
    x_out   = [];
    x_cur   = x0(:);
    frozen  = zeros(4, N_int);
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', dt);

    for j = 1:N_int
        t_j   = t_start + (j-1)*epsilon;
        t_jp1 = min(t_j + epsilon, t_end);
        if t_j >= t_end
            frozen = frozen(:, 1:j-1);
            break
        end

        z_frozen     = x_to_z(x_cur, d0, l1);
        frozen(:, j) = z_frozen;

        odefun = @(t, x) rhs_x(t, x, z_frozen, d0, l1, omega, k12, k112, ...
                               alpha, gamma_val, h_fun);

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
            t_out = t_seg;  x_out = x_seg;
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
        zf    = frozen(:, j_idx);
        u1_out(i) = u1_ctrl(t_out(i), zf, omega, k12, k112, alpha, gamma_val, h_fun);
        u2_out(i) = u2_ctrl(t_out(i), zf, omega, k12, k112, alpha, gamma_val, h_fun);
    end

    % the same trajectory in z, for direct comparison with simulate_case
    z_out = zeros(size(x_out));
    for i = 1:size(x_out, 1)
        z_out(i, :) = x_to_z(x_out(i, :).', d0, l1).';
    end
end

function dx = rhs_x(t, x, zf, d0, l1, omega, k12, k112, alpha, gamma_val, h_fun)
%RHS_X  Exact kinematics (fields_exact), driven by the controller.
    [f1, f2] = fields_exact(x, d0, l1);
    u1 = u1_ctrl(t, zf, omega, k12, k112, alpha, gamma_val, h_fun);
    u2 = u2_ctrl(t, zf, omega, k12, k112, alpha, gamma_val, h_fun);
    dx = f1*u1 + f2*u2;
end

% ---- CONTROLLER -- u1, u2, amplitudes, worst-case bound ----

function dz = rhs_approx(t, z, zf, g1_fun, g2_fun, omega, k12, k112, ...
                         alpha, gamma_val, h_fun)
%RHS_APPROX  Closed-loop right-hand side of the nilpotent approximation,
%   z' = u1*g1(z) + u2*g2(z), with u1, u2 from the controller at the frozen
%   state zf. Used by simulate_case and simulate_case_adaptive.
%
%   g1_fun, g2_fun may be [] -- then the fields come from fields_approx
%   (the normal case). Passing your own handles still works.

    u1 = u1_ctrl(t, zf, omega, k12, k112, alpha, gamma_val, h_fun);
    u2 = u2_ctrl(t, zf, omega, k12, k112, alpha, gamma_val, h_fun);
    if isempty(g1_fun)
        [g1, g2] = fields_approx(z, alpha);
    else
        g1 = g1_fun(z, alpha);
        g2 = g2_fun(z, alpha);
    end
    dz = g1*u1 + g2*u2;
end

function u1_val = u1_ctrl(t, zf, omega, k12, k112, alpha, gamma_val, h_fun)
%U1_CTRL  Driving velocity u1 [m/s] at time t, for the state zf frozen at the
%   start of the current sampling interval. ONE-TRAILER (r=3): two
%   oscillating terms. Amplitudes and term sizes come from osc_terms.
    [a1, ~, ~, ~, T12, T112] = osc_terms(zf, alpha, gamma_val, omega, k12, k112, h_fun);
    u1_val = a1 ...
        + T12  * cos(k12*omega*t) ...
        + T112 * cos(k112*omega*t);
end

function u2_val = u2_ctrl(t, zf, omega, k12, k112, alpha, gamma_val, h_fun)
%U2_CTRL  Steering velocity u2 [rad/s] at time t, for the state zf frozen at
%   the start of the current sampling interval. ONE-TRAILER (r=3): two
%   oscillating terms. Amplitudes and term sizes come from osc_terms.
    [~, a2, s12, s112, T12, T112] = osc_terms(zf, alpha, gamma_val, omega, k12, k112, h_fun);
    u2_val = a2 ...
        + T12 * s12 * sin(k12*omega*t) ...
        - 2*T112 * s112 * cos(2*k112*omega*t);
end

function [a1, a2, s12, s112, T12, T112] = osc_terms( ...
        zf, alpha, gamma_val, omega, k12, k112, h_fun)
%OSC_TERMS  Everything the one-trailer controller needs at the frozen state
%   zf: the constant parts a1, a2, the signs s.. of the bracket amplitudes
%   and the sizes T.. of the two oscillating terms.
%
%   Shared by u1_ctrl / u2_ctrl (the signal) and control_bound (its
%   worst case), so the two can never disagree.

    [a1, a2, a12, a112] = amps(zf, alpha, gamma_val, h_fun);

    T12  = (2*k12*omega*abs(a12))^(1/2);
    T112 = (2*k112*omega)^(2/3) * abs(a112)^(1/3);

    s12 = sign(a12);  s112 = sign(a112);
end

function [a1, a2, a12, a112] = amps(zf, alpha, gamma_val, h_fun)
%AMPS  The four amplitudes, a(z) = -gamma*h(z)*G^{-1}(z)*z. ONE-TRAILER
%   (r=3) version; the closed forms are checked against the symbolic
%   result in derive.m ('symbolic'). Used via osc_terms by the controller and
%   by control_bound.
    h    = h_fun(zf);
    a1   = -gamma_val*h*zf(1);
    a2   = -gamma_val*h*zf(2);
    a12  = -gamma_val*h*(zf(1)*zf(2) + zf(3));
    a112 = -gamma_val*h*(alpha*zf(2)^3 + zf(1)*zf(3) + zf(4));
end

function [b1, b2, terms1, terms2] = control_bound(z, alpha, gamma_val, ...
                                                  omega, k12, k112, h_fun)
%CONTROL_BOUND  Worst-case magnitude of the controls, r=3 (1 trailer).
%   Sums the absolute values of all contributions (peaks-coincide bound,
%   see the 2-trailer control_bound). One fewer term than the 2-trailer
%   version -- no k1112/a1112.
%
%   terms1 = [|a1|, T12, T112]     contributions to |u1|
%   terms2 = [|a2|, T12, 2*T112]   contributions to |u2|

    % h_fun is the state-dependent factor of the Lyapunov function (see
    % params.m). It scales ALL amplitudes, so leaving it out makes this
    % bound wrong by exactly that factor.
    if nargin < 7 || isempty(h_fun)
        h_fun = @(zz) 1;
    end
    [a1, a2, ~, ~, T12, T112] = osc_terms(z, alpha, gamma_val, omega, k12, k112, h_fun);

    terms1 = [abs(a1), T12,   T112];
    terms2 = [abs(a2), T12, 2*T112];

    b1 = sum(terms1);
    b2 = sum(terms2);
end

% ---- GAINS -- admissible gamma and time to converge ----

function g = gamma_max(z, alpha, omega, k12, k112, u1max, u2max, h_fun)
%GAMMA_MAX  Largest gain for which BOTH velocity limits still hold at state z.
%   ONE-TRAILER (r=3) version. Bisection in log space: 200 halvings
%   between 1e-16 and 10.

    if nargin < 8 || isempty(h_fun)
        h_fun = @(zz) 1;
    end

    lo = 1e-16;
    hi = 10;
    for i = 1:200
        mid = sqrt(lo*hi);
        [b1, b2] = control_bound(z, alpha, mid, omega, k12, k112, h_fun);
        if b1 <= u1max && b2 <= u2max
            lo = mid;
        else
            hi = mid;
        end
    end
    g = lo;
end

function [g, hmax] = gamma_ray(z0, alpha, omega, epsilon, k12, k112, ...
                               u1max, u2max, h_fun, cap_product)
%GAMMA_RAY  Largest fixed gain that is admissible along the WHOLE straight
%   ray from z0 to the origin (every radial V drives the state along it).
%
%   [g, hmax] = gamma_ray(P.z0, P.alpha, P.omega, P.epsilon, P.k12, P.k112, ...
%                         P.u1_max, P.u2_max, P.h_fun, P.cap_product)
%
%   Two constraints, both must hold at every point of the ray:
%     (1) amplitude:   control_bound <= (u1max, u2max)      [gamma_max]
%     (2) contraction: epsilon*gamma*h(z) <= cap_product
%   g is the largest gamma meeting both; hmax is the largest h on the ray.
%   This is how the gain in params.m was chosen. Looking at z0 alone is not
%   enough for the sqrt / log V: h grows towards the origin, so (2) can
%   bind near the goal even when (1) is slack at the start.

    r0   = norm(z0);
    zhat = z0(:) / r0;
    rs   = r0 * logspace(-3, 0, 200);

    g = inf;  hmax = 0;
    for r = rs
        g    = min(g, gamma_max(r*zhat, alpha, omega, k12, k112, u1max, u2max, h_fun));
        hmax = max(hmax, h_fun(r*zhat));
    end
    g = min(g, cap_product/(epsilon*hmax));
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

    r0   = norm(z0);
    zhat = z0(:) / r0;
    hr   = @(r) arrayfun(@(rr) h_fun(rr*zhat), r);

    t = zeros(size(fracs));
    for i = 1:numel(fracs)
        t(i) = integral(@(r) 1 ./ (gamma_val*hr(r).*r), fracs(i)*r0, r0);
    end
end

% ---- COORDINATES AND VECTOR FIELDS ----

function [g1, g2] = fields_approx(z, alpha)
%FIELDS_APPROX  Vector fields of the nilpotent approximation of the car with
%   ONE off-hooked trailer (r = 3, 4 states):  z' = u1*g1(z) + u2*g2(z).
%
%   The ONE place these formulas live (calculations, the simulators and
%   the compare_* tasks all use it). derive.m ('symbolic') writes the
%   same fields symbolically.

    g1 = [1; 0; -z(2); -z(3)];
    g2 = [0; 1; 0; -alpha*z(2)^2];
end

function [f1, f2] = fields_exact(x, d0, l1)
%FIELDS_EXACT  Exact kinematics of the car with ONE off-hooked trailer:
%   x' = u1*f1(x) + u2*f2(x),  x = [x1; x2; phi; theta1], angles in radians.
%   These are rows 1..4 of the 5-state two-trailer model (the first four
%   rows never involve the second trailer; see derive.m, 'report').

    x3 = x(3);  x4 = x(4);

    f1 = [cos(x3); sin(x3); 0;  sin(x3 - x4)/l1];
    f2 = [0;       0;       1; -d0*cos(x3 - x4)/l1];
end

function z = x_to_z(x, d0, l1)
%X_TO_Z  Physical state -> privileged coordinates. ONE-TRAILER (r=3).
%
%   x = [x1; x2; phi; theta1]   angles in RADIANS (unlike x0_to_z0 in params.m,
%                               which takes degrees, for convenience when
%                               setting a scenario up by hand)
%
%   This is the version to call inside a loop -- the controller speaks z,
%   the robot measures x, so this runs at the start of every sampling
%   interval. The controls u1, u2 need NO conversion: they are the
%   physical driving / steering velocities and are identical in both
%   descriptions, because x -> z changes only the STATE coordinates.

    x = x(:);
    z = [ x(1);
          x(3);
         -x(2);
          l1*x(2) - d0*l1^2/(d0+l1)*x(3) - l1^3/(d0+l1)*x(4) ];
end

function x = z_to_x(z, d0, l1)
%Z_TO_X  Privileged coordinates -> physical state. Exact inverse of
%   x_to_z (angles in RADIANS). Useful for plotting where the vehicle
%   actually drives, in metres and radians, instead of the abstract z.

    z = z(:);
    A = d0*l1^2/(d0+l1);
    B = l1^3/(d0+l1);
    x = [ z(1);
         -z(3);
          z(2);
          (-l1*z(3) - A*z(2) - z(4))/B ];
end

% ---- k-PAIR SCREEN ----

function ok = heuristic_no_multiple(a, b)
%HEURISTIC_NO_MULTIPLE  Practical placeholder non-resonance check for the
%   2-multiplier (r=3, one-trailer) reduction: neither k12 nor k112 may be
%   an integer multiple of the other.
%
%   THIS IS NOT THE PAPER'S CONDITION C1. That condition is proven only
%   for the 5-state (r=4), 3-multiplier construction -- see the caveat in
%   params.m in this folder. Which condition the 2-multiplier case needs
%   is OPEN (to ask Grushkovskaya); this is just a common-sense stand-in.
%   Do not present results filtered by this as more than a heuristic screen.

    ok = mod(a, b) ~= 0 && mod(b, a) ~= 0;
end
