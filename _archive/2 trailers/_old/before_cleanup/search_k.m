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
                                   U1_MAX, U2_MAX);
                    u1 = control_bound(z0, alpha, kappa, theta, gamma, w, a, b, c);
                    rows(end+1, :) = [a, b, c, g, u1];   %#ok<AGROW>
                end
            end
        end
    end

    [~, idx] = sort(rows(:, 4), 'descend');
    rows = rows(idx, :);

    T = array2table(rows, 'VariableNames', ...
        {'k12', 'k112', 'k1112', 'gamma_max', 'max_u1'});

    base = gamma_max(z0, alpha, kappa, theta, w, 8, 6, 5, U1_MAX, U2_MAX);
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
