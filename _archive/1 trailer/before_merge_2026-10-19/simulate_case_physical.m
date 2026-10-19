function [t_out, x_out, z_out, u1_out, u2_out] = simulate_case_physical( ...
        x0, d0, l1, alpha, gamma_val, epsilon, omega, k12, k112, ...
        t_start, t_end, dt, h_fun)
%SIMULATE_CASE_PHYSICAL  Closed loop on the REAL model. ONE-TRAILER (r=3).
%
%   simulate_case.m integrates the nilpotent approximation
%       z' = u1*g1(z) + u2*g2(z),
%   i.e. the truncated system the controller was designed on. This file
%   integrates the exact kinematics instead,
%       x' = u1*f1(x) + u2*f2(x),
%   driven by exactly the same controller. The gap between the two runs is
%   the remainder the approximation throws away -- which is what the real
%   robot, and Gazebo, will actually feel.
%
%   x0 = [x1; x2; phi; theta1] in RADIANS.
%
%   Loop per sampling interval (also the loop the robot has to run):
%     1. read the physical state x
%     2. z = x_to_z(x)              <- state conversion
%     3. u1, u2 from the controller at that frozen z
%     4. apply u1, u2               <- NO conversion back; u1, u2 are
%                                      already the physical velocities

    % the -1e-9 stops float noise (e.g. 3.0000000000000004) adding an empty interval
    N_int   = ceil((t_end - t_start)/epsilon - 1e-9);
    t_out   = [];
    x_out   = [];
    x_cur   = x0(:);
    frozen  = zeros(4, N_int);
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', dt);

    for j = 1:N_int
        t_j   = t_start + (j-1)*epsilon;
        t_jp1 = min(t_j + epsilon, t_end);
        if t_j >= t_end
            frozen = frozen(:, 1:j-1);
            break
        end

        z_frozen     = x_to_z(x_cur, d0, l1);
        frozen(:, j) = z_frozen;

        odefun = @(t, x) rhs_x(t, x, z_frozen, d0, l1, omega, k12, k112, ...
                               alpha, gamma_val, h_fun);

        % output grid of this interval; it must END exactly at t_jp1,
        % otherwise the next interval starts from a state at an earlier time
        tspan_j = t_j:dt:t_jp1;
        if tspan_j(end) < t_jp1 - 1e-12
            tspan_j = [tspan_j, t_jp1];
        end
        if numel(tspan_j) < 2
            tspan_j = [t_j, t_jp1];
        end

        [t_seg, x_seg] = ode45(odefun, tspan_j, x_cur, options);

        if isempty(t_out)
            t_out = t_seg;  x_out = x_seg;
        else
            t_out = [t_out; t_seg(2:end)];
            x_out = [x_out; x_seg(2:end, :)];
        end
        x_cur = x_seg(end, :).';
    end

    n      = numel(t_out);
    u1_out = zeros(n, 1);
    u2_out = zeros(n, 1);
    n_int  = size(frozen, 2);
    for i = 1:n
        j_idx = min(max(floor((t_out(i) - t_start)/epsilon) + 1, 1), n_int);
        zf    = frozen(:, j_idx);
        u1_out(i) = u1_ctrl(t_out(i), zf, omega, k12, k112, alpha, gamma_val, h_fun);
        u2_out(i) = u2_ctrl(t_out(i), zf, omega, k12, k112, alpha, gamma_val, h_fun);
    end

    % the same trajectory in z, for direct comparison with simulate_case.m
    z_out = zeros(size(x_out));
    for i = 1:size(x_out, 1)
        z_out(i, :) = x_to_z(x_out(i, :).', d0, l1).';
    end
end


function dx = rhs_x(t, x, zf, d0, l1, omega, k12, k112, alpha, gamma_val, h_fun)
%RHS_X  Exact kinematics (fields_exact.m), driven by the controller.
    [f1, f2] = fields_exact(x, d0, l1);
    u1 = u1_ctrl(t, zf, omega, k12, k112, alpha, gamma_val, h_fun);
    u2 = u2_ctrl(t, zf, omega, k12, k112, alpha, gamma_val, h_fun);
    dx = f1*u1 + f2*u2;
end
