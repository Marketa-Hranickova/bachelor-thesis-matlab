
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
%  PART 1 -- symbolic system, Lie brackets, controllability
% =========================================================================
syms z1 z2 z3 z4 z5 real
syms alpha_s kappa_s theta_s real
syms gamma_s positive

z_sym = [z1; z2; z3; z4; z5];

% Vector fields of the nilpotent approximation. The system is
%     z' = u1*g1(z) + u2*g2(z)
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
% Lyapunov function Q(z) = gamma/2 * ||z||^2   (eq. 17)
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
% The controller below uses hand-written closed forms for a(z). Verify that
% they agree with the symbolic result at a random point, so a typo in either
% place is caught immediately.
z_test = [0.37; -0.21; 0.44; -0.13; 0.29];

a_from_sym = double(subs(a_sym, ...
    [z_sym.', gamma_s, alpha_s, kappa_s, theta_s], ...
    [z_test.', gamma, alpha, kappa, theta]));

a_by_hand = [ -gamma*z_test(1); ...
              -gamma*z_test(2); ...
              -gamma*(z_test(1)*z_test(2) + z_test(3)); ...
              -gamma*(z_test(1)*z_test(3) + z_test(4) + alpha*z_test(2)^3); ...
              -gamma*( (alpha + 2*kappa - theta)*z_test(1)*z_test(2)^3 ...
                     + (alpha + 3*kappa)*z_test(2)^2*z_test(3) ...
                     + z_test(1)*z_test(4) + z_test(5) ) ];

res = max(abs(a_from_sym - a_by_hand));
fprintf('\n=== symbolic vs hand-coded a(z): max |difference| = %.3e ===\n', res);
assert(res < 1e-10, 'a(z) mismatch -- symbolic and hand-coded forms disagree');
disp('   PASSED')

% =========================================================================
%  PART 3 -- one closed-loop simulation
% =========================================================================
h_fun  = @(z) 1;
g1_fun = @(z, alpha, kappa, theta) [1; 0; -z(2); -z(3); kappa*z(2)^3 - z(4)];
g2_fun = @(z, alpha, kappa, theta) [0; 1; 0; -alpha*z(2)^2; -theta*z(2)*z(3)];

fprintf('\n========================================\n')
fprintf('  Simulation\n')
fprintf('  z0      = [%.1f; %.1f; %.1f; %.1f; %.1f],  ||z0||^2 = %.2f\n', z0, norm(z0)^2);
fprintf('  time    : [%.1f, %.1f], dt = %.1e\n', t_start, t_end, dt);
fprintf('  epsilon = %.2f, omega = %.4f\n', epsilon, omega);
fprintf('  k       = (%d, %d, %d)\n', k12, k112, k1112);
fprintf('========================================\n')

[t_sim, z_sim, u1_sim, u2_sim] = simulate_case( ...
    z0, alpha, kappa, theta, gamma, ...
    epsilon, omega, k12, k112, k1112, ...
    t_start, t_end, dt, h_fun, g1_fun, g2_fun);

fprintf('  %d steps\n', length(t_sim));
fprintf('  final state : [%.6f; %.6f; %.6f; %.6f; %.6f]\n', z_sim(end,:));
fprintf('  final norm  : %.6f\n', norm(z_sim(end,:)));
fprintf('  max |u1|    : %.3f\n', max(abs(u1_sim)));
fprintf('  max |u2|    : %.3f\n', max(abs(u2_sim)));

% Lyapunov function along the trajectory
V_sim = sum(z_sim.^2, 2);
fprintf('  V(0) = %.6f  ->  V(end) = %.6f\n', V_sim(1), V_sim(end));

fprintf('\n=== Done. ===\n')


% =========================================================================
%  LOCAL FUNCTIONS  (must be last in a MATLAB script)
% =========================================================================
function u1_val = u1_ctrl(t, zf, omega, k12, k112, k1112, ...
                          alpha, kappa, theta, gamma_val, h_fun)
    [a1, ~, a12, a112, a1112] = amps(zf, alpha, kappa, theta, gamma_val, h_fun);
    u1_val = a1 ...
        + (2*k12*omega*abs(a12))^(1/2) * cos(k12*omega*t) ...
        + (2*k112*omega)^(2/3) * abs(a112)^(1/3) * cos(k112*omega*t) ...
        + 2*(k1112*omega)^(3/4) * abs(a1112)^(1/4) * cos(k1112*omega*t);
end


function u2_val = u2_ctrl(t, zf, omega, k12, k112, k1112, ...
                          alpha, kappa, theta, gamma_val, h_fun)
    [~, a2, a12, a112, a1112] = amps(zf, alpha, kappa, theta, gamma_val, h_fun);
    u2_val = a2 ...
        + (2*k12*omega*abs(a12))^(1/2) * sign(a12) * sin(k12*omega*t) ...
        - 2*(2*k112*omega)^(2/3) * abs(a112)^(1/3) * sign(a112) ...
          * cos(2*k112*omega*t) ...
        - 6*(k1112*omega)^(3/4) * abs(a1112)^(1/4) * sign(a1112) ...
          * sin(3*k1112*omega*t);
end


function [a1, a2, a12, a112, a1112] = amps(zf, alpha, kappa, theta, gamma_val, h_fun)
%AMPS  The five amplitudes, a(z) = -gamma*h(z)*G^{-1}(z)*z.
    h     = h_fun(zf);
    a1    = -gamma_val*h*zf(1);
    a2    = -gamma_val*h*zf(2);
    a12   = -gamma_val*h*(zf(1)*zf(2) + zf(3));
    a112  = -gamma_val*h*(zf(1)*zf(3) + zf(4) + alpha*zf(2)^3);
    a1112 = -gamma_val*h*( (alpha + 2*kappa - theta)*zf(1)*zf(2)^3 ...
                         + (alpha + 3*kappa)*zf(2)^2*zf(3) ...
                         + zf(1)*zf(4) + zf(5) );
end


function [t_out, z_out, u1_out, u2_out] = simulate_case( ...
        z0, alpha, kappa, theta, gamma_val, ...
        epsilon, omega, k12, k112, k1112, ...
        t_start, t_end, dt, h_fun, g1_fun, g2_fun)
%SIMULATE_CASE  Sampling scheme: the amplitudes are computed from the state
%   at the start of each interval and held frozen for its whole length.

    N_int   = ceil((t_end - t_start)/epsilon);
    t_out   = [];
    z_out   = [];
    z_cur   = z0(:);
    frozen  = zeros(5, N_int);          % the frozen state of each interval
    t_edges = zeros(1, N_int);
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', dt);

    for j = 1:N_int
        t_j   = t_start + (j-1)*epsilon;
        t_jp1 = min(t_j + epsilon, t_end);
        if t_j >= t_end
            frozen = frozen(:, 1:j-1);
            t_edges = t_edges(1:j-1);
            break
        end

        z_frozen     = z_cur;
        frozen(:, j) = z_frozen;
        t_edges(j)   = t_j;

        odefun = @(t, z) ...
              g1_fun(z, alpha, kappa, theta) * ...
                  u1_ctrl(t, z_frozen, omega, k12, k112, k1112, ...
                          alpha, kappa, theta, gamma_val, h_fun) ...
            + g2_fun(z, alpha, kappa, theta) * ...
                  u2_ctrl(t, z_frozen, omega, k12, k112, k1112, ...
                          alpha, kappa, theta, gamma_val, h_fun);

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

    % Controls along the trajectory. Index straight into the stored frozen
    % states -- the original searched the whole time vector for every sample,
    % which is O(N^2) and unusably slow at dt = 1e-4.
    n      = numel(t_out);
    u1_out = zeros(n, 1);
    u2_out = zeros(n, 1);
    n_int  = size(frozen, 2);
    for i = 1:n
        j_idx = min(max(floor((t_out(i) - t_start)/epsilon) + 1, 1), n_int);
        zf    = frozen(:, j_idx);
        u1_out(i) = u1_ctrl(t_out(i), zf, omega, k12, k112, k1112, ...
                            alpha, kappa, theta, gamma_val, h_fun);
        u2_out(i) = u2_ctrl(t_out(i), zf, omega, k12, k112, k1112, ...
                            alpha, kappa, theta, gamma_val, h_fun);
    end
end
