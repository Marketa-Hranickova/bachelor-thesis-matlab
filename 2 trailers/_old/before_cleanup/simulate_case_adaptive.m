function [t_out, z_out, u1_out, u2_out, gamma_out] = simulate_case_adaptive( ...
        z0, alpha, kappa, theta, ...
        epsilon, omega, k12, k112, k1112, ...
        t_start, t_end, dt, h_fun, g1_fun, g2_fun, ...
        u1_max, u2_max, cap_product)


    N_int     = ceil((t_end - t_start)/epsilon);
    t_out     = [];
    z_out     = [];
    z_cur     = z0(:);
    frozen    = zeros(5, N_int);        % frozen state of each interval
    gamma_int = zeros(1, N_int);        % gain actually used on each interval
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

        % ---- adaptive gain for this interval -------------------------------
        g_hw = gamma_max(z_frozen, alpha, kappa, theta, omega, ...
                          k12, k112, k1112, u1_max, u2_max);
        g_j  = min(g_hw, cap_product/epsilon);
        gamma_int(j) = g_j;

        odefun = @(t, z) ...
              g1_fun(z, alpha, kappa, theta) * ...
                  u1_ctrl(t, z_frozen, omega, k12, k112, k1112, ...
                          alpha, kappa, theta, g_j, h_fun) ...
            + g2_fun(z, alpha, kappa, theta) * ...
                  u2_ctrl(t, z_frozen, omega, k12, k112, k1112, ...
                          alpha, kappa, theta, g_j, h_fun);

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
        u1_out(i) = u1_ctrl(t_out(i), zf, omega, k12, k112, k1112, ...
                            alpha, kappa, theta, g_j, h_fun);
        u2_out(i) = u2_ctrl(t_out(i), zf, omega, k12, k112, k1112, ...
                            alpha, kappa, theta, g_j, h_fun);
    end
end
