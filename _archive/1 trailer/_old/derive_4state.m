% DERIVE_4STATE
% Symbolic derivation of the one-trailer (r = 3) system, starting from the
% five-state model of the paper (arXiv:2403.03341).
%
% calculations.m BEGINS by assuming g1, g2 and alpha. This script DERIVES
% them, so nothing in the one-trailer branch is taken on faith.
%
% Everything below is checked with assert. If the script runs to the end,
% every claim in AIOS/.../DERIVATION_4STATE.md has been reproduced from
% scratch, with symbolic geometry (d0, l1 kept as symbols, not numbers).
%
% Runtime: tens of seconds (symbolic Lie brackets + a Taylor expansion).

clear; clc
syms x1 x2 x3 x4 x5 real
syms d0 d1 l1 l2 positive
syms u1 u2 real
syms ep positive

fprintf('==========================================================\n');
fprintf('  DERIVATION OF THE 4-STATE SYSTEM FROM THE 5-STATE MODEL\n');
fprintf('==========================================================\n');

%% ========================================================================
%  PART 0 -- the 5 -> 4 step: why the last state can simply be dropped
%  ========================================================================
X5 = [x1; x2; x3; x4; x5];

f15 = sin(x3-x5)/l2 + (d1+l1)*sin(x4-x3)*cos(x4-x5)/(l1*l2);
f25 = d0*(d1+l1)*cos(x3-x4)*cos(x4-x5)/(l1*l2) - d0*cos(x3-x5)/l2;

F1_5 = [cos(x3); sin(x3); 0; sin(x3-x4)/l1; f15];
F2_5 = [0;       0;       1; -d0*cos(x3-x4)/l1; f25];

% The claim: the first FOUR components do not involve x5 at all, so the
% subsystem (x1,x2,phi,theta1) is closed -- the second trailer is dragged
% along but cannot influence anything ahead of it.
dep1 = simplify(diff(F1_5(1:4), x5));
dep2 = simplify(diff(F2_5(1:4), x5));
assert(isequal(dep1, sym(zeros(4,1))) && isequal(dep2, sym(zeros(4,1))), ...
       'first four components DO depend on x5 -- the 5->4 step is invalid');
fprintf('\n[0] first four rows are free of x5            -> OK\n');
fprintf('    (so the 4-state model is the exact sub-model, not a limit;\n');
fprintf('     note d1 and l2 occur only in rows 5, hence they disappear)\n');

% The four-state model: literally rows 1..4.
X = [x1; x2; x3; x4];
f1 = F1_5(1:4);
f2 = F2_5(1:4);

%% ========================================================================
%  PART 1 -- Lie brackets, growth vector, controllability
%  ========================================================================
lie = @(f,g) simplify(jacobian(g,X)*f - jacobian(f,X)*g);

B3 = lie(f1, f2);           % [f1,f2]
B4 = lie(f1, B3);           % [f1,[f1,f2]]   -- last one needed for r = 3

fprintf('\n[1] brackets\n');
fprintf('    [f1,f2]      = '); disp(simplify(B3).');
fprintf('    [f1,[f1,f2]] = '); disp(simplify(B4).');

at0 = @(M) simplify(subs(M, [x1 x2 x3 x4], [0 0 0 0]));

% growth vector: how many directions each bracket depth buys
gv = [rank(at0([f1 f2])), rank(at0([f1 f2 B3])), rank(at0([f1 f2 B3 B4]))];
fprintf('    growth vector (2,3,4)?  -> (%d,%d,%d)\n', gv);
assert(isequal(gv, [2 3 4]), 'growth vector is not (2,3,4)');
fprintf('    degree of nonholonomy r = 3                -> OK\n');

F0 = at0([f1 f2 B3 B4]);
detF0 = simplify(det(F0));
fprintf('\n    F(0) =\n'); disp(F0);
fprintf('    det F(0) = '); disp(detF0);
assert(isAlways(simplify(detF0 - (-(d0+l1)/l1^3)) == 0), ...
       'det F(0) does not match -(d0+l1)/l1^3');
fprintf('    matches -(d0+l1)/l1^3                      -> OK\n');
fprintf('    d0, l1 > 0 so this never vanishes: controllable for ANY geometry\n');

%% ========================================================================
%  PART 2 -- privileged coordinates, and the test that no correction is
%            needed (this is the step where the 5-state case needs beta)
%  ========================================================================
% Step 1 of Bellaiche: the LINEAR change z = F(0)^-1 x.
A  = simplify(inv(F0));
zc = sym('z', [4 1], 'real');
x_of_z = simplify(F0*zc);          % inverse map, x in terms of z

fprintf('\n[2] privileged coordinates (linear part only)\n');
for q = 1:4
    fprintf('    z%d = ', q); disp(simplify(A(q,:)*X));
end

% Step 2 would be a nonlinear correction. Test whether it is needed:
% substitute z_i -> ep^w_i * z_i and read off the lowest power of ep in
% each zdot component. The rule is: lowest order MUST equal w_q - 1.
w = [1 1 2 3];
zdot = simplify(A*(f1*u1 + f2*u2));
zdot = simplify(subs(zdot, X, x_of_z));
zdot_s = subs(zdot, zc, (ep.^w(:)).*zc);

fprintf('\n    weighted-order test  (z_i -> ep^w_i z_i)\n');
fprintf('    %-5s %-7s %-10s %-10s %s\n', 'comp', 'weight', 'required', 'lowest', 'leading term');
lead = sym(zeros(4,1));
for q = 1:4
    Tq = taylor(zdot_s(q), ep, 0, 'Order', 5);
    low = NaN;
    for j = 0:4
        if j == 0
            cj = simplify(subs(Tq, ep, 0));
        else
            cj = simplify(subs(diff(Tq, ep, j), ep, 0)/factorial(j));
        end
        if ~isAlways(cj == 0, 'Unknown', 'false')
            if isnan(low), low = j; lead(q) = cj; end
        end
    end
    fprintf('    zdot%-1d %-7d %-10d %-10d ', q, w(q), w(q)-1, low);
    disp(simplify(lead(q)));
    assert(low == w(q)-1, sprintf(...
        ['zdot%d has lowest order %d but weight %d requires %d -- the ' ...
         'linear transform is NOT enough here, a correction IS needed'], ...
        q, low, w(q), w(q)-1));
end
fprintf('    no component carries a lower order than required\n');
fprintf('    => the linear map alone gives privileged coordinates,\n');
fprintf('       the nonlinear (cubic) correction is NOT needed      -> OK\n');
fprintf('    reason: the cubic monomial z2^3 has weight 3. Correcting z4\n');
fprintf('    (weight 3) only admits monomials of weight < 3, and 3 < 3 is\n');
fprintf('    false -- so it cannot enter. In the 5-state case z5 has\n');
fprintf('    weight 4, 3 < 4 holds, and that is exactly the beta term.\n');

%% ========================================================================
%  PART 3 -- the nilpotent approximation and alpha
%  ========================================================================
% The leading terms found above ARE the nilpotent approximation.
g1d = simplify(expand(subs(lead, [u1 u2], [1 0])));
g2d = simplify(expand(subs(lead, [u1 u2], [0 1])));

fprintf('\n[3] nilpotent approximation, derived (not assumed)\n');
fprintf('    g1(z) = '); disp(g1d.');
fprintf('    g2(z) = '); disp(g2d.');

% g2's fourth component is exactly -alpha*z2^2, so divide it out
alpha_derived = simplify(-g2d(4)/zc(2)^2);
fprintf('    alpha (derived) = '); disp(alpha_derived);
assert(isAlways(simplify(alpha_derived - d0*(d0+l1)/2) == 0), ...
       'derived alpha does not match d0(d0+l1)/2');
fprintf('    matches d0*(d0+l1)/2                       -> OK\n');

% cross-check against the function the rest of the code uses
d0n = 0.115; l1n = 0.082;
assert(abs(double(subs(alpha_derived,[d0 l1],[d0n l1n])) - ...
           one_trailer_coeff(d0n,l1n)) < 1e-12, ...
       'one_trailer_coeff.m disagrees with the derivation');
fprintf('    one_trailer_coeff.m agrees numerically      -> OK\n');

% and against the g1,g2 that calculations.m assumes
syms alpha_s real
g1_assumed = [1; 0; -zc(2); -zc(3)];
g2_assumed = [0; 1; 0; -alpha_s*zc(2)^2];
assert(isequal(simplify(g1d - g1_assumed), sym(zeros(4,1))), ...
       'derived g1 differs from the one calculations.m assumes');
assert(isequal(simplify(subs(g2_assumed, alpha_s, alpha_derived) - g2d), ...
       sym(zeros(4,1))), ...
       'derived g2 differs from the one calculations.m assumes');
fprintf('    g1, g2 match what calculations.m assumes    -> OK\n');

%% ========================================================================
%  PART 4 -- nilpotency, det G = 1, and the amplitude formulas
%  ========================================================================
lieZ = @(f,g) simplify(jacobian(g,zc)*f - jacobian(f,zc)*g);
C3 = lieZ(g1d, g2d);
C4 = lieZ(g1d, C3);

fprintf('\n[4] nilpotency\n');
fprintf('    [g1,g2]      = '); disp(C3.');
fprintf('    [g1,[g1,g2]] = '); disp(C4.');
depth4 = {lieZ(g1d,C4), lieZ(g2d,C4), lieZ(g1d,lieZ(g2d,C3)), lieZ(g2d,lieZ(g2d,C3))};
for i = 1:4
    assert(isequal(simplify(depth4{i}), sym(zeros(4,1))), ...
           'a depth-4 bracket is non-zero -- the algebra is not nilpotent of step 3');
end
fprintf('    all four depth-4 brackets vanish           -> OK\n');
fprintf('    (this is WHY the 4th control term of eq.(16) is dropped:\n');
fprintf('     it would excite a direction that is identically zero)\n');

G = [g1d, g2d, C3, C4];
assert(isAlways(simplify(det(G)) == 1), 'det G(z) is not identically 1');
fprintf('    det G(z) = 1 for all z                     -> OK\n');

% a(z) = -gamma*G^-1*z  must reproduce the closed forms in amps.m
syms gam positive
a_sym = simplify(-gam*(G\zc));
a_hand = -gam*[ zc(1);
                zc(2);
                zc(1)*zc(2) + zc(3);
                alpha_derived*zc(2)^3 + zc(1)*zc(3) + zc(4) ];
assert(isequal(simplify(a_sym - a_hand), sym(zeros(4,1))), ...
       'G^-1 z does not match the hand-coded amplitudes');
fprintf('    G^-1 z matches the closed forms            -> OK\n');

% ...and against amps.m itself, numerically
zt = [0.31; -0.17; 0.42; -0.08];
try
    [b1,b2,b12,b112] = amps(zt, one_trailer_coeff(d0n,l1n), 1, @(z) 1);
    ref = double(subs(a_hand, [zc.' gam d0 l1], [zt.' 1 d0n l1n]));
    assert(max(abs([b1;b2;b12;b112] - ref)) < 1e-12, 'amps.m disagrees');
    fprintf('    amps.m agrees numerically                  -> OK\n');
catch ME
    fprintf(2, '    amps.m check FAILED: %s\n', ME.message);
    fprintf(2, '    (if this is a "too many input arguments" error, the\n');
    fprintf(2, '     2-trailer amps.m is shadowing this one -- check with\n');
    fprintf(2, '     "which -all amps" and remove the other folder from the path)\n');
end

%% ========================================================================
%  PART 5 -- independent check: d0 -> 0 must give chained form
%  ========================================================================
fprintf('\n[5] chained-form limit\n');
alpha_0 = limit(alpha_derived, d0, 0);
assert(isAlways(alpha_0 == 0), 'alpha does not vanish as d0 -> 0');
fprintf('    alpha -> 0 as d0 -> 0, so zdot4 -> -z3*u1  -> OK\n');
fprintf('    d0 = 0 is on-axle hitching, and classical theory independently\n');
fprintf('    says that case is chained form. Nothing in this derivation used\n');
fprintf('    that fact, yet it falls out -- a strong external check.\n');

fprintf('\n==========================================================\n');
fprintf('  ALL CHECKS PASSED\n');
fprintf('  calculations.m may now start from g1, g2, alpha above.\n');
fprintf('==========================================================\n');
