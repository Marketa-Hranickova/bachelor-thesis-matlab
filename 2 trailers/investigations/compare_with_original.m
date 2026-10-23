% COMPARE_WITH_ORIGINAL  Runs OUR exact-model simulator (trailer2.m,
% simulate_case_physical) with the settings of the older folder
% original_system (params_o.m) and prints the same quantities as
% calculations_o.m, so the two sets of numbers can be compared line by line.
%
% Usage:
%   1) in the folder original_system run:   calculations_o      (needs the Symbolic Math Toolbox)
%      and keep the lines "Final x: ..." and "||x(end)|| = ..." of Cases 1a, 1b, 2a, 2b
%   2) in this folder (2 trailers/investigations) run:  compare_with_original
%   Same start, same controller, same exact vector fields -> the printed
%   final x must agree (to about 1e-6 or better).

clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));   % trailer2.m and params.m are one folder up
if isappdata(0, 'ctrl_phase'), rmappdata(0, 'ctrl_phase'); end   % no phase shift

gam = 1;  eps_ = 0.1;  om = 2*pi/eps_;
k12 = 8;  k112 = 6;  k1112 = 5;
x0 = [0.4; -0.4; 0.3; -0.2; 0.1];
t0 = 0;  t1 = 5;  dt = 1e-4;

h1 = @(z) 1;
h2 = @(z) abs(log(sum(z.^2)));

% name, d0, d1, l1, l2, alpha, kappa, theta, h
cases = { ...
  '1a', 0.1, 0.1, 1, 1, 0.10495, 0.0024, 0.0112, h1; ...
  '1b', 0.1, 0.1, 1, 1, 0.10495, 0.0024, 0.0112, h2; ...
  '2a', 1,   1,   1, 1, 1,       1.5,    4,      h1; ...
  '2b', 1,   1,   1, 1, 1,       1.5,    4,      h2};

for c = 1:size(cases, 1)
    nm = cases{c,1};  d0 = cases{c,2};  d1 = cases{c,3};
    l1 = cases{c,4};  l2 = cases{c,5};
    al = cases{c,6};  ka = cases{c,7};  th = cases{c,8};  hf = cases{c,9};
    fprintf('\nCase %s (ours) ...\n', nm);
    [~, x, z] = trailer2('simulate_case_physical', x0, d0, d1, l1, l2, ...
        al, ka, th, gam, eps_, om, k12, k112, k1112, t0, t1, dt, hf);
    fprintf('  Final x: [%.6f; %.6f; %.6f; %.6f; %.6f]\n', x(end,:));
    fprintf('  ||x(end)|| = %.6f\n', norm(x(end,:)));
    fprintf('  z(end)  : [%.6f; %.6f; %.6f; %.6f; %.6f]\n', z(end,:));
    V = sum(z(end,:).^2);
    fprintf('  V(end)  = %.6f\n', V);
end
