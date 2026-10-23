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
A = {P.alpha, P.kappa, P.theta};          % the three coefficients, for brevity
K = {P.omega, P.k12, P.k112, P.k1112};    % frequency and multipliers

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
n_fail = n_fail + check('alpha, kappa, theta, beta', ...
    [P.alpha P.kappa P.theta P.beta], ref.coeffs, 1e-12);
n_fail = n_fail + check('z0 = x_to_z(x0)', P.z0(:).', ref.z0, 1e-12);

% ---- 2. controller at z0 -------------------------------------------------
[b1, b2] = control_bound(P.z0, A{:}, P.gamma, K{:}, P.h_fun);
n_fail = n_fail + check('control_bound at z0', [b1 b2], ref.bound, 1e-12);
g = gamma_max(P.z0, A{:}, K{:}, P.u1_max, P.u2_max, P.h_fun);
n_fail = n_fail + check('gamma_max at z0', g, ref.gmax, 1e-10);
[u1, u2] = u_ctrl(0.37, P.z0, K{:}, A{:}, P.gamma, P.h_fun);
n_fail = n_fail + check('u_ctrl(0.37 s) at z0', [u1 u2], ref.uctrl, 1e-12);

% ---- 3. 12 s closed loop, all three simulators ----------------------------
T = {0, 12, P.dt};                         % t_start, t_end, dt
fprintf('\nsimulating 12 s (approx / adaptive / exact) ...\n');

[~, z, u1, u2] = simulate_case(P.z0, A{:}, P.gamma, P.epsilon, K{:}, T{:}, P.h_fun);
n_fail = n_fail + check('approx   z(12 s)', z(end,:), ref.z_approx, 1e-7);
n_fail = n_fail + check('approx   max|u|', [max(abs(u1)) max(abs(u2))], ref.u_approx, 1e-7);

[~, z, u1, u2] = simulate_case_adaptive(P.z0, A{:}, P.epsilon, K{:}, T{:}, P.h_fun, ...
                                        [], [], P.u1_max, P.u2_max, P.cap_product);
n_fail = n_fail + check('adaptive z(12 s)', z(end,:), ref.z_adapt, 1e-7);
n_fail = n_fail + check('adaptive max|u|', [max(abs(u1)) max(abs(u2))], ref.u_adapt, 1e-7);

[~, x, ~, u1, u2] = simulate_case_physical(P.x0, P.d0, P.d1, P.l1, P.l2, A{:}, ...
                                           P.gamma, P.epsilon, K{:}, T{:}, P.h_fun);
n_fail = n_fail + check('exact    x(12 s)', x(end,:), ref.x_exact, 1e-7);
n_fail = n_fail + check('exact    max|u|', [max(abs(u1)) max(abs(u2))], ref.u_exact, 1e-7);

% ---- verdict ---------------------------------------------------------------
fprintf('\n');
if n_fail == 0
    fprintf('=== SELFTEST PASSED ===\n');
else
    fprintf('=== SELFTEST FAILED: %d check(s) moved -- see the FAIL lines above ===\n', n_fail);
end
end


function bad = check(label, got, want, tol)
%CHECK  Compare two numeric vectors; print one line; return 1 if they differ.
    err = max(abs(got(:) - want(:)));
    ok  = err <= tol;
    if ok, tag = 'ok  '; else, tag = 'FAIL'; end
    fprintf('  [%s] %-26s max |difference| = %.2e  (tol %.0e)\n', tag, label, err, tol);
    bad = double(~ok);
end
