% DERIVE_SYMBOLIC  Lie brackets, controllability and the amplitudes a(z) -- symbolically.
%
%   Needs the Symbolic Math Toolbox. Run it once to see (or re-check) the
%   algebra behind the controller; it is NOT needed for simulations.
%
%   Part 1  builds the vector fields g1, g2 of the nilpotent approximation,
%           their Lie brackets up to depth 4, the matrix
%               G(z) = [g1, g2, [g1,g2], [g1,[g1,g2]], [g1,[g1,[g1,g2]]]]
%           and checks that det G = 1 (so G is invertible everywhere).
%   Part 2  derives  a(z) = -G^-1(z) * grad Q(z)  for Q = gamma/2*||z||^2
%           and checks that the closed forms in lib/amps.m agree with it.
%           If you ever change amps.m (or the vector fields), this is the
%           test that catches a typo.
%
%   Only the quadratic V (h = 1) is derived here. Any radial V reduces to
%   the same formulas scaled by h(z) -- see params.m, section 4.

P = params();
alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;  gamma = P.gamma;

fprintf('geometry : d0=%.3f  d1=%.3f  l1=%.3f  l2=%.3f  [m]\n', P.d0, P.d1, P.l1, P.l2);
fprintf('coeffs   : alpha=%.6f  kappa=%.6f  theta=%.6f\n', alpha, kappa, theta);

% =========================================================================
%  PART 1 -- symbolic system, Lie brackets, controllability
% =========================================================================
syms z1 z2 z3 z4 z5 real
syms alpha_s kappa_s theta_s real
syms gamma_s positive

z_sym = [z1; z2; z3; z4; z5];

% Vector fields of the nilpotent approximation. The system is
%     z' = u1*g1(z) + u2*g2(z)
% (the same formulas as lib/fields_approx.m, with symbolic coefficients)
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

% =========================================================================
%  PART 2 -- the amplitudes a(z), for h = 1
% =========================================================================
% Lyapunov function Q(z) = gamma/2 * ||z||^2   (eq. 17 of the paper)
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
% Compare the symbolic result with lib/amps.m (the code every simulation
% actually uses) at a generic point. amps.m is called with the plain
% quadratic V, i.e. h = 1, to match the derivation above.
z_test = [0.37; -0.21; 0.44; -0.13; 0.29];

a_from_sym = double(subs(a_sym, ...
    [z_sym.', gamma_s, alpha_s, kappa_s, theta_s], ...
    [z_test.', gamma, alpha, kappa, theta]));

Pq = P;
Pq.h_fun = lyapunov_h('quadratic');
[b1, b2, b12, b112, b1112] = amps(z_test, Pq, gamma);
a_from_code = [b1; b2; b12; b112; b1112];

res = max(abs(a_from_sym - a_from_code));
fprintf('\n=== symbolic vs lib/amps.m: max |difference| = %.3e ===\n', res);
assert(res < 1e-10, 'a(z) mismatch -- symbolic result and lib/amps.m disagree');
disp('   PASSED')

fprintf('\n=== Done. ===\n')
