function [t_out, z_out, u1_out, u2_out, gamma_out] = simulate_case_adaptive( ...
        z0, alpha, ...
        epsilon, omega, k12, k112, ...
        t_start, t_end, dt, h_fun, g1_fun, g2_fun, ...
        u1_max, u2_max, cap_product)
%SIMULATE_CASE_ADAPTIVE  ONE-TRAILER (r=3) version. Same idea as the
%   2-trailer version: gamma is recomputed at the start of every interval
%   as the largest value that (a) keeps the commanded |u1|,|u2| within the
%   robot's hardware limits at the frozen state (gamma_max.m), and (b)
%   keeps eps*gamma*h(z) under cap_product.
%
%   g1_fun, g2_fun may be [] -- then the fields come from fields_approx.m.
%   gamma_out is the gain in use at every output time.

    % the -1e-9 stops float noise (e.g. 3.0000000000000004) adding an empty interval
    N_int     = ceil((t_end - t_start)/epsilon - 1e-9);
    t_out     = [];
    z_out     = [];
    z_cur     = z0(:);
    frozen    = zeros(4, N_int);
    gamma_int = zeros(1, N_int);
    options   = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', dt);

    for j = 1:N_int
        t_j   = t_start + (j-1)*epsilon;
        t_jp1 = min(t_j + epsilon, t_end);
        if t_j >= t_end
            frozen    = frozen(:, 1:j-1);
            gamma_int = gamma_int(1:j-1);
            break
        end

        z_frozen     = z_cur;
        frozen(:, j) = z_frozen;

        % h_fun must go in here too: the gain is chosen from the same
        % amplitudes the controller then uses, so a mismatch would pick a
        % gain that violates the limits by exactly the factor h.
        g_hw = gamma_max(z_frozen, alpha, omega, k12, k112, u1_max, u2_max, h_fun);
        % contraction cap eps*gamma*h <= cap_product (params.m, rule (2)).
        % For the quadratic V h = 1 and this is the old eps*gamma <= cap;
        % for the sqrt / log V h varies with the state and must be included.
        g_j  = min(g_hw, cap_product/(epsilon*h_fun(z_frozen)));
        gamma_int(j) = g_j;

        odefun = @(t, z) rhs_approx(t, z, z_frozen, g1_fun, g2_fun, omega, ...
                                    k12, k112, alpha, g_j, h_fun);

        % output grid of this interval; it must END exactly at t_jp1,
        % otherwise the next interval starts from a state at an earlier time
        tspan_j = t_j:dt:t_jp1;
        if tspan_j(end) < t_jp1 - 1e-12
            tspan_j = [tspan_j, t_jp1];
        end
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

    n         = numel(t_out);
    u1_out    = zeros(n, 1);
    u2_out    = zeros(n, 1);
    gamma_out = zeros(n, 1);
    n_int     = size(frozen, 2);
    for i = 1:n
        j_idx = min(max(floor((t_out(i) - t_start)/epsilon) + 1, 1), n_int);
        zf    = frozen(:, j_idx);
        g_j   = gamma_int(j_idx);
        gamma_out(i) = g_j;
        u1_out(i) = u1_ctrl(t_out(i), zf, omega, k12, k112, alpha, g_j, h_fun);
        u2_out(i) = u2_ctrl(t_out(i), zf, omega, k12, k112, alpha, g_j, h_fun);
    end
end
