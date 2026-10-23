function R = k_test(opts)
%K_TEST  Does the exact model converge for other multipliers k = (k12,k112,k1112)?
%   bisect_regime step S3 (only k changed, (8,6,5) -> (14,5,1)) was the first
%   step at which the exact model stopped converging. This tests the k's directly.
%
%   k_test                    both parts, three triples each
%   k_test(struct('part','A'))  only the quick part (old regime)
%   k_test(struct('part','B'))  only the design regime (slower)
%   opts.k_list = [8 6 5; 14 5 1; 18 4 1]   rows = triples to test
%   opts.scales = [0.25 0.5 1]   start scales (1 = full design start); default 0.25
%   opts.TmaxB  = 40000           time limit [s] of a part-B run (default 40000)
%
%   Part A  old regime of bisect_regime step S2: robot geometry, start x0.25 (ours), eps = 0.1,
%           gamma = 1, h = |log||z||^2| (no limits). Each run <= 10 s.
%   Part B  design regime: eps = 3 s, sqrt Lyapunov function (delta = 0.005), the SAME start;
%           gamma = the largest the velocity limits and the contraction cap allow FOR THAT k
%           (as in compare_V). Each run stops at ||z|| < 1% (CONVERGED), > 3x (DIVERGED)
%           or after opts.TmaxB (default 40000 s). Slow: minutes per triple.
%   Triples tested must satisfy condition C1 of the paper (no resonance up to order 4);
%   (8,6,5), (14,5,1), (18,4,1) do (checked earlier).
%
%   Columns: min, final = ||z||/||z0||; hitch = largest |angle difference| [rad]; max u1 [m/s].

    if nargin < 1, opts = struct(); end
    if ~isfield(opts, 'part'),   opts.part = 'AB'; end
    if ~isfield(opts, 'k_list'), opts.k_list = [8 6 5; 14 5 1; 18 4 1]; end
    if ~isfield(opts, 'scales'), opts.scales = 0.25; end   % start = scale*(-0.2,-0.2,45,45,45 deg); 1 = full design start
    if ~isfield(opts, 'TmaxB'),  opts.TmaxB = 40000; end   % [s] limit of a part-B run (full start needs ~31000 s)
    setappdata(0, 'ctrl_phase', [0 0 0]);
    addpath(fullfile(fileparts(mfilename('fullpath')), '..'));   % trailer2.m and params.m are one folder up
    Tr = @(name, varargin) trailer2(name, varargin{:});

    geo = struct('d0', 0.115, 'd1', 0.084, 'l1', 0.082, 'l2', 0.089);
    [al, ka, th] = offhooked_coeffs(geo.d0, geo.d1, geo.l1, geo.l2);
    h_log  = @(z) abs(log(sum(z.^2)));
    h_sqrt = @(z) 1/sqrt(norm(z)^2 + 0.005^2);
    U1 = 0.22;  U2 = 2.84;  cap = 0.8;

    R = struct([]);
    for sc = opts.scales(:).'
    x0 = sc*[-0.2; -0.2; deg2rad(45); deg2rad(45); deg2rad(45)];
    z0 = x_to_z(x0, geo.d0, geo.d1, geo.l1, geo.l2);  r0 = norm(z0);
    fprintf('\n########## start scale %g: ||z0|| = %.4f\n', sc, r0);
    for part = opts.part
        if part == 'A', ep = 0.1; hf = h_log;  Tmax = 10;    desc = 'A: eps = 0.1, gamma = 1, h = |log|';
        else,           ep = 3;   hf = h_sqrt; Tmax = opts.TmaxB; desc = 'B: eps = 3, sqrt h, gamma from the limits'; end
        om = 2*pi/ep;
        fprintf('\n===== part %s =====\n', desc);
        fprintf('%-14s %10s %8s %8s %8s %9s %9s %8s  %s\n', 'k', 'gamma', 't_stop', 'min', 'final', 'hitch1', 'hitch2', 'max u1', 'verdict');
        for row = 1:size(opts.k_list, 1)
            k = opts.k_list(row, :);
            if part == 'A'
                gam = 1;
            else
                zhat = z0/r0;  g = inf;  hmax = 0;
                for r = r0*logspace(-3, 0, 200)
                    g    = min(g, Tr('gamma_max', r*zhat, al, ka, th, om, k(1), k(2), k(3), U1, U2, hf));
                    hmax = max(hmax, hf(r*zhat));
                end
                gam = min(g, cap/(ep*hmax));
            end
            dt = min(0.02, 2*pi/(max(k)*om)/40);
            tchunk = ceil(Tmax/40/ep)*ep;
            xc = x0;  t_a = 0;  zmin = 1;  zfin = 1;  h1 = 0;  h2 = 0;  umax = 0;  verdict = 'not finished';
            while t_a < Tmax - 1e-9
                t_b = min(t_a + tchunk, Tmax);
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
            fprintf('(%2d,%2d,%2d)    %10.3e %8.1f %8.4f %8.4f %7.2f r %7.2f r %7.3f   %s\n', ...
                    k, gam, t_a, zmin, zfin, h1, h2, umax, verdict);
            R(end+1).scale = sc;  R(end).part = part;  R(end).k = k;  R(end).gamma = gam;  R(end).t_stop = t_a;
            R(end).zmin = zmin;  R(end).zfin = zfin;  R(end).verdict = verdict;
        end
    end
    end
end
