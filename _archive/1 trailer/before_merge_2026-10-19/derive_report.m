% DERIVE_REPORT
% Runs derive_privileged.m on BOTH systems and checks the output against
% what is independently known.
%
% The engine is told only f1 and f2. It is never told r, the weights, the
% correction, alpha or beta -- it derives them. This script then compares:
%
%   5-state : against the paper (arXiv:2403.03341), whose beta and alpha
%             are published. Reproducing beta from scratch is the real test.
%   4-state : against one_trailer_coeff.m and amps.m, i.e. against the code
%             the rest of this folder runs on.
%
% Runtime: a few minutes, most of it the 5-state case.

clear; clc
syms x1 x2 x3 x4 x5 real

% measured geometry
d0 = sym(115)/1000;  d1 = sym(84)/1000;
l1 = sym(82)/1000;   l2 = sym(89)/1000;

f15 = sin(x3-x5)/l2 + (d1+l1)*sin(x4-x3)*cos(x4-x5)/(l1*l2);
f25 = d0*(d1+l1)*cos(x3-x4)*cos(x4-x5)/(l1*l2) - d0*cos(x3-x5)/l2;
F1  = [cos(x3); sin(x3); 0; sin(x3-x4)/l1; f15];
F2  = [0;       0;       1; -d0*cos(x3-x4)/l1; f25];

%% ========================= 4-STATE =====================================
fprintf('================= 4-STATE (one trailer) =================\n');
R4 = derive_privileged(F1(1:4), F2(1:4), [x1 x2 x3 x4]);
assert(R4.ok, 'weighted-order test failed');
assert(isequal(R4.r, 3) && isequal(R4.growth, [2 3 4]), 'r / growth wrong');
assert(isequal(simplify(R4.corr), sym(zeros(4,1))), ...
       'a correction came out non-zero for the 4-state case');

fprintf('\n  derived g1 = '); disp(R4.g1.');
fprintf('  derived g2 = '); disp(R4.g2.');

alpha4 = simplify(-R4.g2(4)/R4.z(2)^2);
fprintf('  derived alpha = '); disp(alpha4);
assert(isAlways(simplify(alpha4 - d0*(d0+l1)/2) == 0), 'alpha wrong');
fprintf('  == d0*(d0+l1)/2                              -> OK\n');
assert(abs(double(alpha4) - one_trailer_coeff(double(d0),double(l1))) < 1e-12);
fprintf('  one_trailer_coeff.m agrees                   -> OK\n');

% amplitudes against amps.m
zt = [0.31; -0.17; 0.42; -0.08];
ref = double(subs(R4.a, [R4.z.' R4.gam], [zt.' 1]));
try
    [b1,b2,b12,b112] = amps(zt, double(alpha4), 1, @(z) 1);
    assert(max(abs([b1;b2;b12;b112] - ref)) < 1e-12, 'amps.m disagrees');
    fprintf('  amps.m agrees numerically                    -> OK\n');
catch ME
    fprintf(2,'  amps.m check FAILED: %s\n', ME.message);
    fprintf(2,'  (if "too many inputs": the 2-trailer amps.m is shadowing\n');
    fprintf(2,'   this one -- run "which -all amps")\n');
end

%% ========================= 5-STATE =====================================
fprintf('\n================= 5-STATE (two trailers) =================\n');
fprintf('(slower -- this is the one that must reproduce the paper''s beta)\n');
R5 = derive_privileged(F1, F2, [x1 x2 x3 x4 x5]);
assert(R5.ok, 'weighted-order test failed');
assert(isequal(R5.r, 4) && isequal(R5.growth, [2 3 4 5]), 'r / growth wrong');

% only z5 may carry a correction
assert(isequal(simplify(R5.corr(1:4)), sym(zeros(4,1))), ...
       'an unexpected correction appeared on z1..z4');
beta_derived = simplify(R5.corr(5)/R5.z(2)^3);
beta_paper   = d0*l2*(d0^2*(d1+l1)*(l1-l2-d1) + l1^3*(d0+l2))/(6*(d0+l2)*l1^2);

fprintf('\n  derived beta = '); disp(beta_derived);
fprintf('  paper''s beta = '); disp(simplify(beta_paper));
assert(isAlways(simplify(beta_derived - beta_paper) == 0), ...
       'DERIVED BETA DOES NOT MATCH THE PAPER');
fprintf('  exact match                                  -> OK\n');
fprintf('  (z5 = y5 - beta*z2^3, and z2 = x3, so this is eq.(9) of the paper)\n');

alpha5 = simplify(-R5.g2(4)/R5.z(2)^2);
fprintf('\n  derived alpha (5-state) = '); disp(alpha5);
fprintf('  numerically %.10g\n', double(alpha5));

fprintf('\n==========================================================\n');
fprintf('  BOTH SYSTEMS DERIVED FROM f1, f2 ALONE.\n');
fprintf('  The 5-state run reproduced the published beta.\n');
fprintf('==========================================================\n');
