function selftest()
%SELFTEST  ONE-TRAILER folder: does it still compute what it computed before?
%
%   Run after ANY edit to the code (about 20 s). Checks
%     1. alpha and z0, and the x -> z -> x round trip,
%     2. the control bound, gamma_max and one value of u1/u2 at z0,
%     3. a 12 s closed-loop run of all three simulators,
%   against reference values, and prints PASSED or which number moved.
%
%   The references were produced on 18 Oct 2026 with params.m as it was
%   then (gamma = 3.2028e-03, eps = 5, k = (3,2), sqrt V with delta = 0.02,
%   x0 = (-0.2, -0.2, 45, 45 deg)). If you change those ON PURPOSE, this
%   test fails on purpose: check the new numbers, then update the table.

P = params();
n_fail = 0;
K = {P.omega, P.k12, P.k112};

% ---- reference values ----------------------------------------------------
ref.alpha  = 1.132750000000000e-02;
ref.z0     = [-2.000000000000000e-01  7.853981633974483e-01 ...
               2.000000000000000e-01 -2.168101725068444e-02];
ref.bound  = [ 2.119078228329150e-01  3.900679149740466e-01];
ref.gmax   =   3.550395847357964e-03;
ref.uctrl  = [ 1.121128800760134e-01  6.258538533368191e-02];   % at t = 0.37 s

% 12 s runs; the ODE solver may differ in the last digits between MATLAB
% and Octave, hence the looser tolerance for these
ref.z_approx = [-2.474523528968e-01  7.832446761761e-01  2.374795922937e-01 -8.751414946495e-03];
ref.u_approx = [ 2.132600407082e-01  3.902267054966e-01];
ref.z_adapt  = [-2.478831391779e-01  7.804096868508e-01  2.376370100693e-01 -8.390606257002e-03];
ref.u_adapt  = [ 2.200000000000e-01  4.045503847606e-01];
ref.x_exact  = [-2.331515570860e-01 -2.334255561516e-01  7.831352153506e-01 -6.760406833783e-01];
ref.u_exact  = [ 2.119078228329e-01  3.900679149740e-01];

% ---- 1. geometry-only numbers and the coordinate map ----------------------
n_fail = n_fail + check('alpha', P.alpha, ref.alpha, 1e-12);
n_fail = n_fail + check('z0 = x0_to_z0(x0)', P.z0(:).', ref.z0, 1e-12);

xr = [P.x0(1); P.x0(2); deg2rad(P.x0(3)); deg2rad(P.x0(4))];
n_fail = n_fail + check('x -> z -> x round trip', ...
    z_to_x(x_to_z(xr, P.d0, P.l1), P.d0, P.l1).', xr.', 1e-12);
n_fail = n_fail + check('x_to_z = x0_to_z0', x_to_z(xr, P.d0, P.l1).', P.z0(:).', 1e-12);

% ---- 2. controller at z0 ---------------------------------------------------
[b1, b2] = control_bound(P.z0, P.alpha, P.gamma, K{:}, P.h_fun);
n_fail = n_fail + check('control_bound at z0', [b1 b2], ref.bound, 1e-12);
g = gamma_max(P.z0, P.alpha, K{:}, P.u1_max, P.u2_max, P.h_fun);
n_fail = n_fail + check('gamma_max at z0', g, ref.gmax, 1e-10);
u1 = u1_ctrl(0.37, P.z0, K{:}, P.alpha, P.gamma, P.h_fun);
u2 = u2_ctrl(0.37, P.z0, K{:}, P.alpha, P.gamma, P.h_fun);
n_fail = n_fail + check('u1/u2 at t = 0.37 s', [u1 u2], ref.uctrl, 1e-12);

% ---- 3. 12 s closed loop, all three simulators -------------------------------
T = {0, 12, P.dt};                         % t_start, t_end, dt
fprintf('\nsimulating 12 s (approx / adaptive / exact) ...\n');

[~, z, u1, u2] = simulate_case(P.z0, P.alpha, P.gamma, P.epsilon, K{:}, T{:}, P.h_fun);
n_fail = n_fail + check('approx   z(12 s)', z(end,:), ref.z_approx, 1e-7);
n_fail = n_fail + check('approx   max|u|', [max(abs(u1)) max(abs(u2))], ref.u_approx, 1e-7);

[~, z, u1, u2] = simulate_case_adaptive(P.z0, P.alpha, P.epsilon, K{:}, T{:}, P.h_fun, ...
                                        [], [], P.u1_max, P.u2_max, P.cap_product);
n_fail = n_fail + check('adaptive z(12 s)', z(end,:), ref.z_adapt, 1e-7);
n_fail = n_fail + check('adaptive max|u|', [max(abs(u1)) max(abs(u2))], ref.u_adapt, 1e-7);

[~, x, ~, u1, u2] = simulate_case_physical(xr, P.d0, P.l1, P.alpha, P.gamma, ...
                                           P.epsilon, K{:}, T{:}, P.h_fun);
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
