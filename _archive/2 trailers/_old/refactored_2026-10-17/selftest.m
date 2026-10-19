function selftest()
%SELFTEST  Does this folder still compute what it computed before?
%
%   Run this after ANY edit to the code (about 30 s). It checks
%     1. the geometry-only numbers (alpha, kappa, theta, beta, z0),
%     2. the coordinate transform (det F(0) = eq. (4) of the paper, and
%        x -> z -> x round trip),
%     3. the control bound and gamma_max at the start state,
%     4. a short 12 s closed-loop run of all three variants (approximate /
%        adaptive / exact) against stored reference values.
%   and prints PASSED or says exactly which number moved.
%
%   The reference values were produced with the parameters in params.m as
%   of 17 Oct 2026 (geometry measured 15 Sep, gamma = 2.7498e-05, eps = 3,
%   k = (14,5,1), sqrt-V with delta = 0.005, x0 = (-0.2,-0.2,45,45,45 deg)).
%   If you change those ON PURPOSE, this test will fail on purpose: run it
%   once, check the new numbers make sense, and update the table below.

P = params();
n_fail = 0;

% ---- reference values (from the pre-refactor code, 17 Oct 2026) ---------
ref.coeffs = [7.355682212287460e-03  2.329616385637499e-02 ...
              6.111464836543101e-02 -1.085630195308095e-04];
ref.z0     = [-2.000000000000000e-01  7.853981633974517e-01 ...
               1.999999999999998e-01 -5.143399189943022e-02  2.492339780897482e-03];
ref.bound  = [ 2.200003688803968e-01  5.335988881155652e-01];
ref.gmax   =   2.749784700008027e-05;

% 12 s run (4 intervals). The ODE solver may differ in the last digits
% between MATLAB and Octave, hence the looser tolerance for these.
ref.z_approx = [-1.999210939499e-01  7.850883000634e-01  1.999172538071e-01 ...
                -5.141390261616e-02  2.485627563645e-03];
ref.u_approx = [ 2.200003688804e-01  5.247694189120e-01];
ref.z_adapt  = [-1.999210855864e-01  7.850882672199e-01  1.999172449357e-01 ...
                -5.141390049781e-02  2.485626801677e-03];
ref.u_adapt  = [ 2.200000000000e-01  5.247685880535e-01];
ref.x_exact  = [-1.999595205010e-01 -1.999301954083e-01  7.850883015238e-01 ...
                 8.395618896470e-01  4.836880264283e-01];
ref.u_exact  = [ 2.200003688804e-01  5.247694189120e-01];

% ---- 1. geometry-only numbers -------------------------------------------
n_fail = n_fail + check('alpha, kappa, theta, beta', ...
    [P.alpha P.kappa P.theta P.beta], ref.coeffs, 1e-12);
n_fail = n_fail + check('z0 = x_to_z(x0)', P.z0(:).', ref.z0, 1e-12);

% ---- 2. coordinate transform --------------------------------------------
check_transform(P);

% ---- 3. control bound and gamma_max at z0 --------------------------------
[b1, b2] = control_bound(P.z0, P, P.gamma);
n_fail = n_fail + check('control_bound at z0', [b1 b2], ref.bound, 1e-12);
n_fail = n_fail + check('gamma_max at z0', gamma_max(P.z0, P), ref.gmax, 1e-10);

% ---- 4. a short closed-loop run, all three variants ----------------------
Q = P;
Q.t_end = 12;

fprintf('\nsimulating 12 s (approx / adaptive / exact) ...\n');
A = simulate(Q, 'approx', 'fixed');
n_fail = n_fail + check('approx   z(12 s)', A.z(end, :), ref.z_approx, 1e-7);
n_fail = n_fail + check('approx   max|u|', [max(abs(A.u1)) max(abs(A.u2))], ref.u_approx, 1e-7);

D = simulate(Q, 'approx', 'adaptive');
n_fail = n_fail + check('adaptive z(12 s)', D.z(end, :), ref.z_adapt, 1e-7);
n_fail = n_fail + check('adaptive max|u|', [max(abs(D.u1)) max(abs(D.u2))], ref.u_adapt, 1e-7);

E = simulate(Q, 'exact', 'fixed');
n_fail = n_fail + check('exact    x(12 s)', E.x(end, :), ref.x_exact, 1e-7);
n_fail = n_fail + check('exact    max|u|', [max(abs(E.u1)) max(abs(E.u2))], ref.u_exact, 1e-7);

% ---- verdict -------------------------------------------------------------
fprintf('\n');
if n_fail == 0
    fprintf('=== SELFTEST PASSED ===\n');
else
    fprintf('=== SELFTEST FAILED: %d check(s) moved -- see the lines marked FAIL above ===\n', n_fail);
end
end


function bad = check(label, got, want, tol)
%CHECK  Compare two numeric vectors; print one line; return 1 if they differ.
    err = max(abs(got(:) - want(:)));
    ok  = err <= tol;
    if ok, tag = 'ok  '; else, tag = 'FAIL'; end
    fprintf('  [%s] %-24s max |difference| = %.2e  (tol %.0e)\n', tag, label, err, tol);
    bad = double(~ok);
end
