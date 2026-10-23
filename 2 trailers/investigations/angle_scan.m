function R = angle_scan(opts)
%ANGLE_SCAN  Does the INITIAL ANGLE matter for convergence of the exact model?
%   Same design as k_test part B (two trailers, eps = 3 s, sqrt Lyapunov
%   function with delta = 0.005, gamma = the largest the velocity limits and
%   the contraction cap allow), same exact kinematics; only the start angles change.
%
%   angle_scan                      default scan
%   angle_scan(struct('pos',[-0.05;-0.05]))   other start position [m]
%   R = angle_scan(...)             one struct per run
%
%   opts.pos        start position (x1; x2) [m]                     default (-0.10; -0.10)
%   opts.angles     angles tested [deg]                              default [0 15 30 45 60]
%   opts.families   1 = phi = theta1 = theta2 = a (no initial hitch angle; the
%                       design start is a = 45 deg, which points at the goal)
%                   2 = phi = a, theta1 = theta2 = 0 (trailers lie along the x axis,
%                       initial hitch angle a)                       default [1 2]
%   opts.k_list     rows = (k12,k112,k1112)                          default [14 5 1; 18 4 1]
%   opts.Tmax       time limit of one run [s]                        default 20000
%
%   A run stops at ||z|| < 1% of ||z0|| (CONVERGED), > 3x (DIVERGED), or at Tmax.
%   Columns: min, final = ||z||/||z0||; hitch1 = car-trailer1, hitch2 = trailer1-trailer2:
%   largest |angle difference| over the run [rad]; max u1 [m/s] (limit 0.22).
%   The URDF joint limit is 1.40 rad (placeholder); the real joint turns about 1.57 rad (90 deg) or less.

    if nargin < 1, opts = struct(); end
    if ~isfield(opts, 'pos'),      opts.pos = [-0.10; -0.10]; end
    if ~isfield(opts, 'angles'),   opts.angles = [0 15 30 45 60]; end
    if ~isfield(opts, 'families'), opts.families = [1 2]; end
    if ~isfield(opts, 'k_list'),   opts.k_list = [14 5 1; 18 4 1]; end
    if ~isfield(opts, 'Tmax'),     opts.Tmax = 20000; end
    setappdata(0, 'ctrl_phase', [0 0 0]);
    addpath(fullfile(fileparts(mfilename('fullpath')), '..'));   % trailer2.m and params.m are one folder up
    Tr = @(name, varargin) trailer2(name, varargin{:});

    geo = struct('d0', 0.115, 'd1', 0.084, 'l1', 0.082, 'l2', 0.089);
    [al, ka, th] = offhooked_coeffs(geo.d0, geo.d1, geo.l1, geo.l2);
    hf = @(z) 1/sqrt(norm(z)^2 + 0.005^2);
    ep = 3;  om = 2*pi/ep;  U1 = 0.22;  U2 = 2.84;  cap = 0.8;
    fname = {'phi=th1=th2=a', 'phi=a, th1=th2=0'};

    fprintf('start position (%.3f, %.3f) m, eps = %g s, sqrt h, Tmax = %g s\n', opts.pos, ep, opts.Tmax);
    fprintf('%-17s %5s %-10s %10s %8s %8s %8s %8s %8s %8s  %s\n', 'family', 'a', 'k', 'gamma', ...
            't_stop', '|z0|', 'min', 'final', 'hitch1', 'hitch2', 'verdict (max u1)');
    R = struct([]);
    for k = opts.k_list.'
        for fam = opts.families
            for a = opts.angles
                ar = deg2rad(a);
                if fam == 1, x0 = [opts.pos(:); ar; ar; ar];
                else,        x0 = [opts.pos(:); ar; 0; 0]; end
                z0 = x_to_z(x0, geo.d0, geo.d1, geo.l1, geo.l2);  r0 = norm(z0);

                zhat = z0/r0;  g = inf;  hmax = 0;               % gamma as in compare_V / k_test
                for r = r0*logspace(-3, 0, 200)
                    g    = min(g, Tr('gamma_max', r*zhat, al, ka, th, om, k(1), k(2), k(3), U1, U2, hf));
                    hmax = max(hmax, hf(r*zhat));
                end
                gam = min(g, cap/(ep*hmax));

                dt = min(0.02, 2*pi/(max(k)*om)/40);
                tchunk = ceil(opts.Tmax/40/ep)*ep;
                xc = x0;  t_a = 0;  zmin = 1;  zfin = 1;  h1 = 0;  h2 = 0;  umax = 0;  verdict = 'not finished';
                while t_a < opts.Tmax - 1e-9
                    t_b = min(t_a + tchunk, opts.Tmax);
                    [~, xs, zs, u1] = Tr('simulate_case_physical', xc, geo.d0, geo.d1, geo.l1, geo.l2, al, ka, th, ...
                                         gam, ep, om, k(1), k(2), k(3), t_a, t_b, dt, hf);
                    r    = sqrt(sum(zs.^2, 2))/r0;
                    zmin = min(zmin, min(r));  zfin = r(end);
                    h1 = max(h1, max(abs(xs(:,3) - xs(:,4))));
                    h2 = max(h2, max(abs(xs(:,4) - xs(:,5))));
                    umax = max(umax, max(abs(u1)));
                    xc = xs(end, :).';  t_a = t_b;
                    if zfin < 0.01, verdict = 'CONVERGED'; break, end
                    if zfin > 3,    verdict = 'DIVERGED';  break, end
                end
                fprintf('%-17s %5g (%2d,%d,%d) %10.3e %8.0f %8.4f %8.4f %8.4f %6.2f r %6.2f r  %s (%.3f)\n', ...
                        fname{fam}, a, k, gam, t_a, r0, zmin, zfin, h1, h2, verdict, umax);
                R(end+1).family = fam;  R(end).angle = a;  R(end).k = k.';  R(end).gamma = gam;
                R(end).t_stop = t_a;  R(end).zmin = zmin;  R(end).zfin = zfin;
                R(end).hitch = [h1 h2];  R(end).umax = umax;  R(end).verdict = verdict;
            end
        end
    end
end
