function out = simulate(P, model, gain_mode)
%SIMULATE  Closed-loop simulation with the sampled oscillating controller.
%
%   out = simulate(P)                        approximate (z) model, fixed gain
%   out = simulate(P, 'exact')               exact kinematics x' = u1 f1 + u2 f2
%   out = simulate(P, 'approx', 'adaptive')  gain re-chosen every interval
%
%   model       'approx' (default)  integrate the nilpotent approximation
%                                   z' = u1*g1(z) + u2*g2(z), start at P.z0
%                                   -- the system the controller is designed on
%               'exact'             integrate the real kinematics, start at
%                                   P.x0 (radians)
%                                   -- what the robot / Gazebo actually is
%   gain_mode   'fixed' (default)   gamma = P.gamma for the whole run
%               'adaptive'          every interval, gamma = the largest value
%                                   that respects the velocity limits at the
%                                   frozen state (gamma_max.m), capped so
%                                   that epsilon*gamma <= P.cap_product
%
%   Both models are driven by EXACTLY the same controller and the same
%   u1, u2 -- only the plant differs. The loop is the one the robot runs:
%       1. read the state (for 'exact': z = x_to_z(x) -- state conversion)
%       2. freeze the amplitudes at that state          (ctrl_frame.m)
%       3. apply u1(t), u2(t) for one interval epsilon  (u_ctrl.m)
%          -- NO conversion of u back: u1, u2 are already the physical
%          driving and steering velocities.
%
%   Result (struct):
%       out.t        n-by-1   time [s]
%       out.z        n-by-5   state in privileged coordinates
%       out.x        n-by-5   physical state (only for 'exact', else [])
%       out.u1,u2    n-by-1   control inputs
%       out.gamma    n-by-1   gain in use at each sample
%       out.model, out.gain_mode   labels

    if nargin < 2 || isempty(model),     model     = 'approx'; end
    if nargin < 3 || isempty(gain_mode), gain_mode = 'fixed';  end

    switch model
        case 'approx', exact = false;
        case 'exact',  exact = true;
        otherwise, error('simulate:model', 'model must be ''approx'' or ''exact''.');
    end
    if ~any(strcmp(gain_mode, {'fixed', 'adaptive'}))
        error('simulate:gain', 'gain_mode must be ''fixed'' or ''adaptive''.');
    end

    eps_s = P.epsilon;
    % number of intervals; the 1e-9 keeps a ratio like 3.0000000000000004
    % from creating a spurious extra interval of length 1e-16
    N_int = ceil((P.t_end - P.t_start)/eps_s - 1e-9);

    if exact, s_cur = P.x0(:); else, s_cur = P.z0(:); end

    t_out     = [];
    s_out     = [];
    frames    = cell(1, N_int);         % frozen controller data of each interval
    gamma_int = zeros(1, N_int);        % gain used on each interval
    options   = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', P.dt);

    % ---- integrate interval by interval ----------------------------------
    for j = 1:N_int
        t_j   = P.t_start + (j-1)*eps_s;
        t_jp1 = min(t_j + eps_s, P.t_end);
        if t_j >= P.t_end
            frames    = frames(1:j-1);
            gamma_int = gamma_int(1:j-1);
            break
        end

        % the controller speaks z, the plant speaks x (exact) or z (approx)
        if exact
            z_frozen = x_to_z(s_cur, P.d0, P.d1, P.l1, P.l2);
        else
            z_frozen = s_cur;
        end

        if strcmp(gain_mode, 'adaptive')
            % gamma_max must see the same h(z) the controller uses, or it
            % would pick a gain that breaks the limits by exactly that factor
            g_j = min(gamma_max(z_frozen, P), P.cap_product/eps_s);
        else
            g_j = P.gamma;
        end
        gamma_int(j) = g_j;

        % amplitudes are constant on the interval: compute them once
        c = ctrl_frame(z_frozen, P, g_j);
        frames{j} = c;

        odefun = @(t, s) rhs(t, s, c, P, exact);

        tspan_j = t_j:P.dt:t_jp1;
        if tspan_j(end) < t_jp1
            tspan_j(end+1) = t_jp1;     %#ok<AGROW> always end ON the boundary
        end
        if numel(tspan_j) < 2
            tspan_j = [t_j, t_jp1];
        end

        [t_seg, s_seg] = ode45(odefun, tspan_j, s_cur, options);

        if isempty(t_out)
            t_out = t_seg;
            s_out = s_seg;
        else
            t_out = [t_out; t_seg(2:end)];          %#ok<AGROW>
            s_out = [s_out; s_seg(2:end, :)];       %#ok<AGROW>
        end
        s_cur = s_seg(end, :).';
    end

    % ---- reconstruct the controls on the output grid ---------------------
    n      = numel(t_out);
    u1     = zeros(n, 1);
    u2     = zeros(n, 1);
    gam    = zeros(n, 1);
    n_int  = numel(frames);
    for i = 1:n
        j_idx  = min(max(floor((t_out(i) - P.t_start)/eps_s) + 1, 1), n_int);
        gam(i) = gamma_int(j_idx);
        [u1(i), u2(i)] = u_ctrl(t_out(i), frames{j_idx}, P);
    end

    % ---- pack the result --------------------------------------------------
    out.t     = t_out;
    out.u1    = u1;
    out.u2    = u2;
    out.gamma = gam;
    out.model = model;
    out.gain_mode = gain_mode;
    if exact
        out.x = s_out;
        out.z = x_to_z(s_out.', P.d0, P.d1, P.l1, P.l2).';
    else
        out.x = [];
        out.z = s_out;
    end
end


function ds = rhs(t, s, c, P, exact)
%RHS  Right-hand side  s' = u1*f1(s) + u2*f2(s)  of the chosen plant.
    [u1, u2] = u_ctrl(t, c, P);
    if exact
        [f1, f2] = fields_exact(s, P);
    else
        [f1, f2] = fields_approx(s, P);
    end
    ds = f1*u1 + f2*u2;
end
