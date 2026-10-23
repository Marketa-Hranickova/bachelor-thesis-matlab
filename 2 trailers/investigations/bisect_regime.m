function R = bisect_regime(steps)
%BISECT_REGIME  Walk from the settings of the OLD project (original_system,
%   where the exact model converged) to OUR design, changing ONE block per
%   step, and report where the exact kinematics stops converging.
%
%   bisect_regime            all steps 0..6
%   bisect_regime(2:4)       only some steps
%   R = bisect_regime(...)   table of results (one struct per step)
%
%   Steps (each keeps everything of the step before):
%     0  old project: paper geometry (26), eps = 0.1, gamma = 1, k = (8,6,5),
%        h = |log||z||^2|, start x0 = (0.4,-0.4,0.3,-0.2,0.1)        [= case 1b]
%     1  + robot geometry; the start POSITIONS scaled by l1 = 0.082 (angles unchanged),
%        so the start is as far from the goal, in units of the trailer length, as before
%     2  + our start, scaled by 0.25:  0.25*(-0.2,-0.2,45,45,45 deg)
%     3  + our multipliers k = (14,5,1)
%     4  + our Lyapunov function, sqrt with delta = 0.005; gamma lowered to 0.04 = 0.8/(eps*max h), because
%        h reaches 1/delta = 200 and eps*gamma*h must stay below ~0.8 (with gamma = 1 it would be 20,
%        i.e. unstable by construction, whatever the geometry or k)
%     5  + eps = 3 s, gamma = 1.3e-3 (largest the CONTRACTION cap allows; velocity limits ignored)
%     6  + gamma = 4.6833e-5 (velocity limits respected) = our design from the x0.25 start
%
%   Each run stops as soon as ||z|| < 1% of ||z0|| (CONVERGED), or > 3x ||z0||
%   (DIVERGED), or at Tmax. Velocity limits are NOT enforced in steps 0-5, the
%   column "max u1" shows how far over the limit 0.22 m/s the run goes.
%   Uses the exact kinematics (trailer2 / simulate_case_physical) only.
%
%   Cost: steps 0-4 are short, steps 5-6 take about a minute each.

    if nargin < 1, steps = 0:6; end
    setappdata(0, 'ctrl_phase', [0 0 0]);               % no phase shift
    addpath(fullfile(fileparts(mfilename('fullpath')), '..'));   % trailer2.m and params.m are one folder up
    Tr = @(name, varargin) trailer2(name, varargin{:});

    paper = struct('d0', 0.1,   'd1', 0.1,   'l1', 1,     'l2', 1);
    robot = struct('d0', 0.115, 'd1', 0.084, 'l1', 0.082, 'l2', 0.089);
    coef_paper = [0.10495 0.0024 0.0112];               % as printed in the paper / params_o.m

    h_log  = @(z) abs(log(sum(z.^2)));
    h_sqrt = @(z) 1/sqrt(norm(z)^2 + 0.005^2);

    x0_old    = [0.4; -0.4; 0.3; -0.2; 0.1];
    x0_scaled = [0.082*x0_old(1:2); x0_old(3:5)];
    x0_ours   = 0.25*[-0.2; -0.2; deg2rad(45); deg2rad(45); deg2rad(45)];

    % name, geometry, coefficients ('paper' or 'robot'), x0, eps, gamma, k, h, Tmax
    S = { ...
      'S0 old project (case 1b)',                paper, 'paper', x0_old,    0.1, 1,       [8 6 5],   h_log,  10;   ...
      'S1 + robot geometry, start scaled',       robot, 'robot', x0_scaled, 0.1, 1,       [8 6 5],   h_log,  10;   ...
      'S2 + our start (x0.25)',                  robot, 'robot', x0_ours,   0.1, 1,       [8 6 5],   h_log,  10;   ...
      'S3 + k = (14,5,1)',                       robot, 'robot', x0_ours,   0.1, 1,       [14 5 1],  h_log,  10;   ...
      'S4 + sqrt Lyapunov function, gamma = 0.04', robot, 'robot', x0_ours, 0.1, 0.04,   [14 5 1],  h_sqrt, 60;   ...
      'S5 + eps = 3 s, gamma = 1.3e-3',          robot, 'robot', x0_ours,   3,   1.3e-3,  [14 5 1],  h_sqrt, 600;  ...
      'S6 + gamma = 4.68e-5 (limits) = ours',    robot, 'robot', x0_ours,   3,   4.6833e-5,[14 5 1], h_sqrt, 6000};

    fprintf('%-40s %8s %8s %8s %8s %9s %9s %8s  %s\n', 'step', '|z0|', 't_stop', 'min', 'final', ...
            'hitch1', 'hitch2', 'max u1', 'verdict');
    R = struct([]);
    for i = steps(:).'
        c = S(i+1, :);
        geo = c{2};  x0 = c{4};  ep = c{5};  gam = c{6};  k = c{7};  hf = c{8};  Tmax = c{9};
        if strcmp(c{3}, 'paper')
            al = coef_paper(1);  ka = coef_paper(2);  th = coef_paper(3);
        else
            [al, ka, th] = offhooked_coeffs(geo.d0, geo.d1, geo.l1, geo.l2);
        end
        om = 2*pi/ep;
        dt = min(0.02, 2*pi/(max(k)*om)/40);
        tchunk = ceil(Tmax/40/ep)*ep;                   % a whole number of intervals
        z0n = norm(x_to_z(x0, geo.d0, geo.d1, geo.l1, geo.l2));

        x = x0(:);  t_a = 0;  zmin = 1;  zfin = 1;  h1 = 0;  h2 = 0;  umax = 0;
        verdict = 'not finished';
        while t_a < Tmax - 1e-9
            t_b = min(t_a + tchunk, Tmax);
            [~, xs, zs, u1, ~] = Tr('simulate_case_physical', x, geo.d0, geo.d1, geo.l1, geo.l2, ...
                                    al, ka, th, gam, ep, om, k(1), k(2), k(3), t_a, t_b, dt, hf);
            r    = sqrt(sum(zs.^2, 2))/z0n;
            zmin = min(zmin, min(r));  zfin = r(end);
            h1 = max(h1, max(abs(xs(:,3) - xs(:,4))));
            h2 = max(h2, max(abs(xs(:,4) - xs(:,5))));
            umax = max(umax, max(abs(u1)));
            x = xs(end, :).';  t_a = t_b;
            if zfin < 0.01, verdict = 'CONVERGED'; break, end
            if zfin > 3,    verdict = 'DIVERGED';  break, end
        end
        fprintf('%-40s %8.4f %8.1f %8.4f %8.4f %7.2f r %7.2f r %7.3f   %s\n', c{1}, z0n, t_a, ...
                zmin, zfin, h1, h2, umax, verdict);
        R(end+1).step = i;  R(end).name = c{1};  R(end).t_stop = t_a;  R(end).zmin = zmin;
        R(end).zfin = zfin;  R(end).hitch = [h1 h2];  R(end).umax = umax;  R(end).verdict = verdict;
    end
    fprintf('\nmin, final = ||z||/||z0||; hitch = largest |angle difference| over all samples [rad]; max u1 in m/s (limit 0.22).\n');
end
