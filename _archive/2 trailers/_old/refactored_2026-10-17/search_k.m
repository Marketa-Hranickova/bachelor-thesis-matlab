function T = search_k(kmax, eps_override)
%SEARCH_K  All C1-admissible multiplier triples, ranked by admissible gain.
%
%   T = search_k(15)        search k = 1..15, epsilon from params.m
%   T = search_k(15, 8.0)   with epsilon = 8 s
%
%   For every triple (k12, k112, k1112) that satisfies condition C1 of the
%   theorem (satisfies_c1.m) it computes
%       gamma_max : the largest gain that keeps |u1|, |u2| inside the robot's
%                   limits at the start state P.z0
%       max_u1    : the worst-case |u1| of that triple at the CURRENT P.gamma
%   and returns a table sorted by gamma_max (best first); the best ten are
%   printed next to the paper's (8,6,5).
%
%   What usually pays is a SMALL k1112: the deepest bracket term carries
%   omega^(3/4) and a factor 6 in u2, so it dominates the demand.
%
%   Slow for large kmax (it loops over kmax^3 triples); 15 is enough.

    if nargin < 1, kmax = 30; end

    P = params();
    if nargin >= 2
        P = set_epsilon(P, eps_override);       % epsilon AND omega
    end
    print_setup(P);

    rows = zeros(0, 5);
    for a = 1:kmax
        for b = 1:kmax
            for c = 1:kmax
                if satisfies_c1(a, b, c)
                    Q = P;
                    Q.k12 = a;  Q.k112 = b;  Q.k1112 = c;
                    g  = gamma_max(P.z0, Q);
                    u1 = control_bound(P.z0, Q, P.gamma);
                    rows(end+1, :) = [a, b, c, g, u1];   %#ok<AGROW>
                end
            end
        end
    end

    [~, idx] = sort(rows(:, 4), 'descend');
    rows = rows(idx, :);

    T = array2table(rows, 'VariableNames', ...
        {'k12', 'k112', 'k1112', 'gamma_max', 'max_u1'});

    % the paper's choice, as the yardstick
    Qp = P;
    Qp.k12 = 8;  Qp.k112 = 6;  Qp.k1112 = 5;
    base = gamma_max(P.z0, Qp);

    fprintf('\n %d admissible triples with k <= %d,  eps = %.2f s\n', ...
            size(rows, 1), kmax, P.epsilon);
    fprintf(' the paper''s (8,6,5) gives gamma_max = %.3e\n\n', base);
    fprintf('%6s %6s %8s %14s %12s %10s\n', ...
            'k12', 'k112', 'k1112', 'gamma_max', 'max|u1|', 'vs paper');
    for i = 1:min(10, size(rows, 1))
        fprintf('%6d %6d %8d %14.4e %12.3f %9.2fx\n', ...
                rows(i,1), rows(i,2), rows(i,3), rows(i,4), rows(i,5), ...
                rows(i,4)/base);
    end
    fprintf('\n');
end
