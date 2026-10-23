function [t_out, z_out, u1_out, u2_out] = simulate_case( ...
        z0, alpha, kappa, theta, gamma_val, ...
        epsilon, omega, k12, k112, k1112, ...
        t_start, t_end, dt, h_fun, g1_fun, g2_fun)
%SIMULATE_CASE  Closed loop on the NILPOTENT APPROXIMATION
%   z' = u1*g1(z) + u2*g2(z), fixed gain gamma_val.
%
%   Every epsilon seconds the state is frozen, the controller amplitudes
%   are computed from it (u_ctrl.m), and ode45 integrates one interval.
%
%   g1_fun, g2_fun are optional: leave them out (or pass []) and the fields
%   come from fields_approx.m.
%   For the exact kinematics use simulate_case_physical.m, for the
%   state-dependent gain simulate_case_adaptive.m.

    if nargin < 15, g1_fun = []; end
    if nargin < 16, g2_fun = []; end

    % the -1e-9 stops float noise (e.g. 3.0000000000000004) adding an empty interval
    N_int   = ceil((t_end - t_start)/epsilon - 1e-9);
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

        odefun = @(t, z) rhs_approx(t, z, z_frozen, g1_fun, g2_fun, omega, ...
                               k12, k112, k1112, alpha, kappa, theta, ...
                               gamma_val, h_fun);

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

    n      = numel(t_out);
    u1_out = zeros(n, 1);
    u2_out = zeros(n, 1);
    n_int  = size(frozen, 2);
    for i = 1:n
        j_idx = min(max(floor((t_out(i) - t_start)/epsilon) + 1, 1), n_int);
        zf    = frozen(:, j_idx);
        [u1_out(i), u2_out(i)] = u_ctrl(t_out(i), zf, omega, k12, k112, ...
                                        k1112, alpha, kappa, theta, gamma_val, h_fun);
    end
end

