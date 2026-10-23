% CALCULATIONS  One closed-loop simulation with the numbers from params.m.
%
%   Runs the approximate (nilpotent) model at the fixed gain and leaves
%   t_sim, z_sim, u1_sim, u2_sim, normz_sim in the workspace for plots.m.
%   No Symbolic Toolbox needed: the symbolic derivation that used to be
%   PART 1-2 of this file is now derive_symbolic.m.

P = params();
d0 = P.d0;  d1 = P.d1;  l1 = P.l1;  l2 = P.l2;
alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;  beta = P.beta;
gamma = P.gamma;  epsilon = P.epsilon;  omega = P.omega;
k12 = P.k12;  k112 = P.k112;  k1112 = P.k1112;
z0 = P.z0;  t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;
U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

fprintf('geometry : d0=%.3f  d1=%.3f  l1=%.3f  l2=%.3f  [m]\n', d0, d1, l1, l2);
fprintf('coeffs   : alpha=%.6f  kappa=%.6f  theta=%.6f\n', alpha, kappa, theta);

% =========================================================================
%  One closed-loop simulation (approximate model, fixed gain)
% =========================================================================
h_fun  = P.h_fun;    % Lyapunov function lives in params.m

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

% state norm along the trajectory (every radial V is a function of it)
normz_sim = sqrt(sum(z_sim.^2, 2));
fprintf('  ||z(0)|| = %.6f  ->  ||z(end)|| = %.6f\n', normz_sim(1), normz_sim(end));

fprintf('\n=== Done. ===\n')
