function T = search_k(kmax, eps_override)
%SEARCH_K  All C1-admissible multiplier triples, ranked by admissible gain.
%
%   T = search_k(15)        search k = 1..15, epsilon from params.m
%   T = search_k(15, 8.0)   with epsilon = 8 s
%
%   Returns a table; prints the best ten. What usually pays is a SMALL
%   k1112: the deepest bracket term carries omega^(3/4) and a factor 6
%   in u2, so it dominates the demand.

    if nargin < 1, kmax = 30; end
    P = params();
    d0 = P.d0;  d1 = P.d1;  l1 = P.l1;  l2 = P.l2;
    alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;  beta = P.beta;
    gamma = P.gamma;  epsilon = P.epsilon;  omega = P.omega;
    % the loop below supplies its own multipliers
    z0 = P.z0;  t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;
    U1_MAX = P.u1_max;  U2_MAX = P.u2_max;
    h_fun  = P.h_fun;

    fprintf('geometry : d0=%.3f  d1=%.3f  l1=%.3f  l2=%.3f  [m]\n', d0, d1, l1, l2);
    fprintf('coeffs   : alpha=%.6f  kappa=%.6f  theta=%.6f\n', alpha, kappa, theta);
    if nargin >= 2
        epsilon = eps_override;
    end
    w = 2*pi/epsilon;


    rows = zeros(0, 5);
    for a = 1:kmax
        for b = 1:kmax
            for c = 1:kmax
                if satisfies_c1(a, b, c)
                    g  = gamma_max(z0, alpha, kappa, theta, w, a, b, c, ...
                                   U1_MAX, U2_MAX, h_fun);
                    u1 = control_bound(z0, alpha, kappa, theta, gamma, w, ...
                                       a, b, c, h_fun);
                    rows(end+1, :) = [a, b, c, g, u1];   %#ok<AGROW>
                end
            end
        end
    end

    [~, idx] = sort(rows(:, 4), 'descend');
    rows = rows(idx, :);

    T = array2table(rows, 'VariableNames', ...
        {'k12', 'k112', 'k1112', 'gamma_max', 'max_u1'});

    base = gamma_max(z0, alpha, kappa, theta, w, 8, 6, 5, ...
                     U1_MAX, U2_MAX, h_fun);
    fprintf('\n %d admissible triples with k <= %d,  eps = %.2f s\n', ...
            size(rows, 1), kmax, epsilon);
    fprintf(' the paper''s (8,6,5) gives gamma_max = %.3e\n\n', base);
    fprintf('%6s %6s %8s %14s %12s %10s\n', ...
            'k12', 'k112', 'k1112', 'gamma_max', 'max|u1|', 'vs paper');
    for i = 1:min(10, size(rows, 1))
        fprintf('%6d %6d %8d %14.4e %12.1f %9.2fx\n', ...
                rows(i,1), rows(i,2), rows(i,3), rows(i,4), rows(i,5), ...
                rows(i,4)/base);
    end
    fprintf('\n');
end


function ok = satisfies_c1(k12, k112, k1112)
%SATISFIES_C1  Condition C1 of Theorem 1 (Zuyev & Grushkovskaya, 2403.03341).
%
%   C1) there are no resonances of order up to 4 in each of the tuples
%       (k12, k112, k1112), (k12, 2*k112, k1112),
%       (k12, k112, 3*k1112), (k12, 2*k112, 3*k1112)
%
%   A tuple (k_1,...,k_s) is in resonance of order K if there exist SETWISE
%   COPRIME integers (c_1,...,c_s) with sum|c_i| = K and sum c_i*k_i = 0.
%
%   Smaller multipliers mean smaller commanded velocities, so it is tempting
%   to shrink them -- this is the condition that says how far you may go.

    tuples = [k12, k112,   k1112;
              k12, 2*k112, k1112;
              k12, k112,   3*k1112;
              k12, 2*k112, 3*k1112];

    ok = true;
    for r = 1:size(tuples, 1)
        if has_resonance(tuples(r, :), 4)
            ok = false;
            return
        end
    end
end


function res = has_resonance(kt, max_order)
%HAS_RESONANCE  True if the tuple has a resonance of order <= max_order.
    res = false;
    for c1 = -max_order:max_order
        for c2 = -max_order:max_order
            for c3 = -max_order:max_order
                K = abs(c1) + abs(c2) + abs(c3);
                if K == 0 || K > max_order
                    continue
                end
                if gcd(gcd(abs(c1), abs(c2)), abs(c3)) ~= 1
                    continue        % not setwise coprime
                end
                if c1*kt(1) + c2*kt(2) + c3*kt(3) == 0
                    res = true;
                    return
                end
            end
        end
    end
end
