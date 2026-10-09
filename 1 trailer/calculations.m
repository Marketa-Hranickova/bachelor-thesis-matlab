% calculations.m -- symbolic derivation and one closed-loop simulation
% of the off-hooked ONE-TRAILER system (r = 3, 4 states).
%
%   Derivation: AIOS/.../DERIVATION_4STATE.md (same method as the paper --
%   Bellaiche/Jean/Risler privileged coordinates -- redone from scratch for
%   this 4-state system, since the paper itself only treats r = 4).
%
% One geometry, one gain function (h = 1), one simulation.
% Local functions are at the END of this file, as MATLAB requires.

P = params();
d0 = P.d0;  l1 = P.l1;
alpha = P.alpha;
gamma = P.gamma;  epsilon = P.epsilon;  omega = P.omega;
k12 = P.k12;  k112 = P.k112;
z0 = P.z0;  t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;
U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

fprintf('geometry : d0=%.3f  l1=%.3f  [m]\n', d0, l1);
fprintf('coeff    : alpha=%.6f\n', alpha);

% =========================================================================
%  PART 1 -- symbolic system, Lie brackets, controllability
% =========================================================================
syms z1 z2 z3 z4 real
syms alpha_s real
syms gamma_s positive

z_sym = [z1; z2; z3; z4];

% Vector fields of the nilpotent approximation (r = 3). The system is
%     z' = u1*g1(z) + u2*g2(z)
% One state fewer than the 2-trailer case -- no z5, no kappa, no theta.
g1_sym = [1; 0; -z2; -z3];
g2_sym = [0; 1; 0; -alpha_s*z2^2];

fprintf('\n=== Vector field g1(z) ===\n'); disp(g1_sym)
fprintf('\n=== Vector field g2(z) ===\n'); disp(g2_sym)

dg1dz = jacobian(g1_sym, z_sym);
dg2dz = jacobian(g2_sym, z_sym);

% [g1, g2]
Lie_g1g2_sym = simplify(dg2dz*g1_sym - dg1dz*g2_sym);
fprintf('\n=== [g1, g2] ===\n'); disp(Lie_g1g2_sym)

% [g1, [g1, g2]] -- this is the LAST bracket needed for r = 3. There is no
% [g1,[g1,[g1,g2]]] here (that belongs to a 5th state that doesn't exist).
dLie_g1g2_dz    = jacobian(Lie_g1g2_sym, z_sym);
Lie_g1_g1g2_sym = simplify(dLie_g1g2_dz*g1_sym - dg1dz*Lie_g1g2_sym);
fprintf('\n=== [g1, [g1, g2]] ===\n'); disp(Lie_g1_g1g2_sym)

% The four brackets span the tangent space -> the system is controllable.
G_sym = [g1_sym, g2_sym, Lie_g1g2_sym, Lie_g1_g1g2_sym];
fprintf('\n=== G matrix ===\n'); disp(G_sym)

detG = simplify(det(G_sym));
fprintf('\n=== det(G) ===\n'); disp(detG)
fprintf('(DERIVATION_4STATE.md says this should be exactly 1, for all z)\n');

% =========================================================================
%  PART 2 -- the amplitudes a(z), for h = 1
% =========================================================================
% Lyapunov function Q(z) = gamma/2 * ||z||^2, now over 4 states.
Q_sym     = gamma_s/2 * (z1^2 + z2^2 + z3^2 + z4^2);
gradQ_sym = jacobian(Q_sym, z_sym).';
fprintf('\n=== grad Q(z) ===\n'); disp(gradQ_sym)

G_inv_sym = simplify(inv(G_sym));

% a(z) = -h * G^{-1} * grad Q,  with h = 1
a_sym = simplify(-G_inv_sym*gradQ_sym);

fprintf('\n=== a(z) components, h = 1 ===\n')
labels = {'a1   ', 'a2   ', 'a12  ', 'a112 '};
for i = 1:4
    fprintf('%s = ', labels{i}); disp(a_sym(i))
end

% ------------------------------------------------------------------ check
% The controller below uses hand-written closed forms for a(z). Verify that
% they agree with the symbolic result at a random point, so a typo in either
% place is caught immediately.
z_test = [0.37; -0.21; 0.44; -0.13];

a_from_sym = double(subs(a_sym, ...
    [z_sym.', gamma_s, alpha_s], ...
    [z_test.', gamma, alpha]));

a_by_hand = [ -gamma*z_test(1); ...
              -gamma*z_test(2); ...
              -gamma*(z_test(1)*z_test(2) + z_test(3)); ...
              -gamma*(alpha*z_test(2)^3 + z_test(1)*z_test(3) + z_test(4)) ];

res = max(abs(a_from_sym - a_by_hand));
fprintf('\n=== symbolic vs hand-coded a(z): max |difference| = %.3e ===\n', res);
assert(res < 1e-10, 'a(z) mismatch -- symbolic and hand-coded forms disagree');
disp('   PASSED')

% =========================================================================
%  PART 3 -- one closed-loop simulation
% =========================================================================
% CAUTION: the control law below (u1_ctrl/u2_ctrl at the end of this file)
% is the first three levels of eq. (16) from the paper (a1, a2, a12, a112,
% same roots/frequencies as the 2-trailer case), with the fourth level
% (a1112, which needs a 5th state) simply dropped. This is the natural
% guess, not something the paper proves for r = 3 -- Theorem 1 is proven
% for the specific 5-state (r = 4) construction. Whether the SAME roots,
% frequencies and condition C1 are actually valid/sufficient for this
% 4-state system needs confirming with the supervisor before the numbers
% below go into the thesis as more than a preliminary estimate.
h_fun  = @(z) 1;
g1_fun = @(z, alpha) [1; 0; -z(2); -z(3)];
g2_fun = @(z, alpha) [0; 1; 0; -alpha*z(2)^2];

fprintf('\n========================================\n')
fprintf('  Simulation\n')
fprintf('  z0      = [%.1f; %.1f; %.1f; %.1f],  ||z0||^2 = %.2f\n', z0, norm(z0)^2);
fprintf('  time    : [%.1f, %.1f], dt = %.1e\n', t_start, t_end, dt);
fprintf('  epsilon = %.2f, omega = %.4f\n', epsilon, omega);
fprintf('  k       = (%d, %d)   [k12, k112 -- no k1112 for r=3]\n', k12, k112);
fprintf('========================================\n')

[t_sim, z_sim, u1_sim, u2_sim] = simulate_case( ...
    z0, alpha, gamma, ...
    epsilon, omega, k12, k112, ...
    t_start, t_end, dt, h_fun, g1_fun, g2_fun);

fprintf('  %d steps\n', length(t_sim));
fprintf('  final state : [%.6f; %.6f; %.6f; %.6f]\n', z_sim(end,:));
fprintf('  final norm  : %.6f\n', norm(z_sim(end,:)));
fprintf('  max |u1|    : %.3f  (limit %.2f)\n', max(abs(u1_sim)), U1_MAX);
fprintf('  max |u2|    : %.3f  (limit %.2f)\n', max(abs(u2_sim)), U2_MAX);

% Lyapunov function along the trajectory
V_sim = sum(z_sim.^2, 2);
fprintf('  V(0) = %.6f  ->  V(end) = %.6f\n', V_sim(1), V_sim(end));

fprintf('\n=== Done. ===\n')


% =========================================================================
%  LOCAL FUNCTIONS  (must be last in a MATLAB script)
% =========================================================================
function u1_val = u1_ctrl(t, zf, omega, k12, k112, alpha, gamma_val, h_fun)
    [a1, ~, a12, a112] = amps(zf, alpha, gamma_val, h_fun);
    u1_val = a1 ...
        + (2*k12*omega*abs(a12))^(1/2) * cos(k12*omega*t) ...
        + (2*k112*omega)^(2/3) * abs(a112)^(1/3) * cos(k112*omega*t);
end


function u2_val = u2_ctrl(t, zf, omega, k12, k112, alpha, gamma_val, h_fun)
    [~, a2, a12, a112] = amps(zf, alpha, gamma_val, h_fun);
    u2_val = a2 ...
        + (2*k12*omega*abs(a12))^(1/2) * sign(a12) * sin(k12*omega*t) ...
        - 2*(2*k112*omega)^(2/3) * abs(a112)^(1/3) * sign(a112) ...
          * cos(2*k112*omega*t);
end


function [a1, a2, a12, a112] = amps(zf, alpha, gamma_val, h_fun)
%AMPS  The four amplitudes, a(z) = -gamma*h(z)*G^{-1}(z)*z, for r = 3.
    h    = h_fun(zf);
    a1   = -gamma_val*h*zf(1);
    a2   = -gamma_val*h*zf(2);
    a12  = -gamma_val*h*(zf(1)*zf(2) + zf(3));
    a112 = -gamma_val*h*(alpha*zf(2)^3 + zf(1)*zf(3) + zf(4));
end


function [t_out, z_out, u1_out, u2_out] = simulate_case( ...
        z0, alpha, gamma_val, ...
        epsilon, omega, k12, k112, ...
        t_start, t_end, dt, h_fun, g1_fun, g2_fun)
%SIMULATE_CASE  Sampling scheme: the amplitudes are computed from the state
%   at the start of each interval and held frozen for its whole length.
%   Same scheme as the 2-trailer version, just 4 states and 2 multipliers.

    N_int   = ceil((t_end - t_start)/epsilon);
    t_out   = [];
    z_out   = [];
    z_cur   = z0(:);
    frozen  = zeros(4, N_int);          % the frozen state of each interval
    t_edges = zeros(1, N_int);
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', dt);

    for j = 1:N_int
        t_j   = t_start + (j-1)*epsilon;
        t_jp1 = min(t_j + epsilon, t_end);
        if t_j >= t_end
            frozen  = frozen(:, 1:j-1);
            t_edges = t_edges(1:j-1);
            break
        end

        z_frozen     = z_cur;
        frozen(:, j) = z_frozen;
        t_edges(j)   = t_j;

        odefun = @(t, z) ...
              g1_fun(z, alpha) * ...
                  u1_ctrl(t, z_frozen, omega, k12, k112, alpha, gamma_val, h_fun) ...
            + g2_fun(z, alpha) * ...
                  u2_ctrl(t, z_frozen, omega, k12, k112, alpha, gamma_val, h_fun);

        tspan_j = t_j:dt:t_jp1;
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

    % Controls along the trajectory -- index straight into the stored frozen
    % states (see the 2-trailer version for why: searching the whole time
    % vector per sample is O(N^2) and unusably slow at dt = 1e-4).
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
