function [t_out, x_out, z_out, u1_out, u2_out] = simulate_case_physical( ...
        x0, d0, d1, l1, l2, alpha, kappa, theta, gamma_val, ...
        epsilon, omega, k12, k112, k1112, ...
        t_start, t_end, dt, h_fun)
%SIMULATE_CASE_PHYSICAL  Closed loop on the REAL model.
%
%   simulate_case.m integrates the nilpotent approximation
%       z' = u1*g1(z) + u2*g2(z),
%   i.e. the truncated system the controller was designed on. This file
%   integrates the exact kinematics instead,
%       x' = u1*f1(x) + u2*f2(x),
%   while driving it with exactly the same controller. The difference
%   between the two runs is the remainder r(z,u) that the approximation
%   throws away -- which is what the real robot will actually feel.
%
%   Loop structure per sampling interval (and this is also the loop the
%   robot has to run):
%       1. read the physical state x
%       2. z = x_to_z(x)                 <- state conversion
%       3. u1, u2 from the controller at that frozen z
%       4. apply u1, u2 to the plant     <- NO conversion back to x:
%                                           u1, u2 are already the
%                                           physical driving / steering
%                                           velocities
%
%   Returns the trajectory in both coordinate sets, plus the controls.

    N_int   = ceil((t_end - t_start)/epsilon);
    t_out   = [];
    x_out   = [];
    x_cur   = x0(:);
    frozen  = zeros(5, N_int);          % frozen z of each interval
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', dt);

    for j = 1:N_int
        t_j   = t_start + (j-1)*epsilon;
        t_jp1 = min(t_j + epsilon, t_end);
        if t_j >= t_end
            frozen = frozen(:, 1:j-1);
            break
        end

        % The only new ingredient: the controller speaks z, the plant
        % speaks x, so convert the measured state once per interval.
        z_frozen     = x_to_z(x_cur, d0, d1, l1, l2);
        frozen(:, j) = z_frozen;

        odefun = @(t, x) phys_rhs(t, x, z_frozen, d0, d1, l1, l2, ...
                                  omega, k12, k112, k1112, ...
                                  alpha, kappa, theta, gamma_val, h_fun);

        tspan_j = t_j:dt:t_jp1;
        if numel(tspan_j) < 2
            tspan_j = [t_j, t_jp1];
        end

        [t_seg, x_seg] = ode45(odefun, tspan_j, x_cur, options);

        if isempty(t_out)
            t_out = t_seg;
            x_out = x_seg;
        else
            t_out = [t_out; t_seg(2:end)];
            x_out = [x_out; x_seg(2:end, :)];
        end
        x_cur = x_seg(end, :).';
    end

    % Controls along the trajectory, indexed straight into the frozen
    % states (same trick as simulate_case.m).
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

    % The same trajectory expressed in z (vectorised form of x_to_z), so
    % it can be compared directly with simulate_case.m.
    [~, ~, ~, beta] = offhooked_coeffs(d0, d1, l1, l2);
    M0         = privileged_frame(d0, d1, l1, l2);
    Y          = (M0 \ x_out.').';
    z_out      = Y;
    z_out(:,5) = Y(:,5) - beta*x_out(:,3).^3;
end


function dx = phys_rhs(t, x, zf, d0, d1, l1, l2, ...
                       omega, k12, k112, k1112, ...
                       alpha, kappa, theta, gamma_val, h_fun)
    [f1, f2] = phys_fields(x, d0, d1, l1, l2);
    u1 = u1_ctrl(t, zf, omega, k12, k112, k1112, ...
                 alpha, kappa, theta, gamma_val, h_fun);
    u2 = u2_ctrl(t, zf, omega, k12, k112, k1112, ...
                 alpha, kappa, theta, gamma_val, h_fun);
    dx = f1*u1 + f2*u2;
end
