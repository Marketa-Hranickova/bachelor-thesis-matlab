function try_k(k12, k112, k1112, eps_override)
%TRY_K  Is this triple of frequency multipliers allowed, and is it any good?
%
%   try_k(7, 4, 2)          uses the epsilon from params.m
%   try_k(7, 4, 2, 8.0)     overrides it with epsilon = 8 s
%
%   Prints the C1 verdict, the worst-case control demand and the largest
%   admissible gain, alongside the paper's own choice (8, 6, 5).

    P = params();
    d0 = P.d0;  d1 = P.d1;  l1 = P.l1;  l2 = P.l2;
    alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;  beta = P.beta;
    gamma = P.gamma;  epsilon = P.epsilon;  omega = P.omega;
    % k12, k112, k1112 come from the ARGUMENTS -- do not unpack them here
    z0 = P.z0;  t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;
    U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

    fprintf('geometry : d0=%.3f  d1=%.3f  l1=%.3f  l2=%.3f  [m]\n', d0, d1, l1, l2);
    fprintf('coeffs   : alpha=%.6f  kappa=%.6f  theta=%.6f\n', alpha, kappa, theta);
    if nargin >= 4
        epsilon = eps_override;
    end
    w = 2*pi/epsilon;


    ok = satisfies_c1(k12, k112, k1112);

    [b1, b2] = control_bound(z0, alpha, kappa, theta, gamma, w, k12, k112, k1112);
    gm = gamma_max(z0, alpha, kappa, theta, w, k12, k112, k1112, U1_MAX, U2_MAX);

    [rb1, rb2] = control_bound(z0, alpha, kappa, theta, gamma, w, 8, 6, 5);
    rgm = gamma_max(z0, alpha, kappa, theta, w, 8, 6, 5, U1_MAX, U2_MAX);

    fprintf('\n k = (%d, %d, %d)   at eps = %.2f s\n', k12, k112, k1112, epsilon);
    fprintf(' ----------------------------------------------\n');
    if ok
        fprintf('  condition C1 : SATISFIED -- this triple may be used\n');
    else
        fprintf('  condition C1 : VIOLATED  -- NOT usable, whatever the numbers below say\n');
    end
    fprintf('  max |u1|     : %10.2f m/s    (limit %.2f  -> %.0f x over)\n', ...
            b1, U1_MAX, b1/U1_MAX);
    fprintf('  max |u2|     : %10.2f rad/s  (limit %.2f  -> %.0f x over)\n', ...
            b2, U2_MAX, b2/U2_MAX);
    fprintf('  gamma_max    : %10.3e\n', gm);
    fprintf('\n  the paper''s (8, 6, 5) at the same eps:\n');
    fprintf('  max |u1|     : %10.2f      gamma_max : %.3e\n', rb1, rgm);
    fprintf('  -> this triple is %.2fx %s on gain\n\n', ...
            max(gm/rgm, rgm/gm), ternary(gm >= rgm, 'BETTER', 'worse'));
end


function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
