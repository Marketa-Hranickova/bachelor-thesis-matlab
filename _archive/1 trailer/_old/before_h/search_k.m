function T = search_k(kmax, eps_override)
%SEARCH_K  Heuristically-admissible (k12,k112) pairs, ranked by admissible
%   gain. ONE-TRAILER (r=3) version -- mirrors "2 trailers/search_k.m",
%   but with only two multipliers and a HEURISTIC screen instead of the
%   paper's proven condition C1 (see heuristic_no_multiple.m: this is NOT
%   confirmed to be the right condition for this reduced case).
%
%   T = search_k(15)        search k = 1..15, epsilon from params.m
%   T = search_k(15, 8.0)   with epsilon = 8 s

    if nargin < 1, kmax = 15; end
    P = params();
    d0 = P.d0;  l1 = P.l1;
    alpha = P.alpha;
    gamma = P.gamma;  epsilon = P.epsilon;  omega = P.omega;
    z0 = P.z0;
    U1_MAX = P.u1_max;  U2_MAX = P.u2_max;

    fprintf('geometry : d0=%.3f  l1=%.3f  [m]\n', d0, l1);
    fprintf('coeff    : alpha=%.6f\n', alpha);
    if nargin >= 2
        epsilon = eps_override;
    end
    w = 2*pi/epsilon;

    rows = zeros(0, 4);
    for a = 1:kmax
        for b = 1:kmax
            if heuristic_no_multiple(a, b)
                g  = gamma_max(z0, alpha, w, a, b, U1_MAX, U2_MAX);
                u1 = control_bound(z0, alpha, gamma, w, a, b);
                rows(end+1, :) = [a, b, g, u1];   %#ok<AGROW>
            end
        end
    end

    [~, idx] = sort(rows(:, 3), 'descend');
    rows = rows(idx, :);

    T = array2table(rows, 'VariableNames', ...
        {'k12', 'k112', 'gamma_max', 'max_u1'});

    base = gamma_max(z0, alpha, w, P.k12, P.k112, U1_MAX, U2_MAX);
    fprintf('\n %d heuristic-admissible pairs with k <= %d,  eps = %.2f s\n', ...
            size(rows, 1), kmax, epsilon);
    fprintf(' current default (%d,%d) gives gamma_max = %.3e\n\n', ...
            P.k12, P.k112, base);
    fprintf('%6s %6s %14s %12s %10s\n', ...
            'k12', 'k112', 'gamma_max', 'max|u1|', 'vs default');
    for i = 1:min(10, size(rows, 1))
        fprintf('%6d %6d %14.4e %12.1f %9.2fx\n', ...
                rows(i,1), rows(i,2), rows(i,3), rows(i,4), rows(i,3)/base);
    end
    fprintf('\n(condition used is a heuristic, not the proven C1 -- see\n');
    fprintf('heuristic_no_multiple.m and the caveat in calculations.m)\n\n');
end
