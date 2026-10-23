function R = paper_real(opts)
%PAPER_REAL  The "real" (exact-kinematics) version of the paper's simulation:
%   parameters of Zuyev & Grushkovskaya, ignoring the robot's own numbers.
%
%   paper_real                 all 8 cases (about 10-20 min in MATLAB)
%   paper_real(struct('T', 3)) shorter runs
%   paper_real(struct('cases', [1 2]))  only some cases
%   R = paper_real(...)        one struct per case
%
%   Common to all cases: eps = 0.1 s, gamma = 1, k = (8,6,5) (paper), no
%   limits on u1, u2. Cases = geometry x start x Lyapunov factor h:
%     geometry   (26): d0 = d1 = 0.1, l1 = l2 = 1   (alpha, kappa, theta as printed in the paper)
%                (27): all four = 1
%     start      "paper": z0 = (1,-2,2,-1,1), the start of the paper's figures, converted
%                         to x0 by z_to_x (exact inverse of x_to_z)
%                "small": x0 = (0.4,-0.4,0.3,-0.2,0.1), the start of the older project original_system
%     h          1 (quadratic V) or |log||z||^2| (the older project's second choice)
%
%   Per case the exact model and the nilpotent approximation are driven by the
%   same controller; ||z||/||z0|| is printed every second for both. A case
%   stops early when the exact model has ||z|| < 1% of ||z0|| (CONVERGED) or
%   > 3x ||z0|| (DIVERGED), else at opts.T (default 8 s).
%
%   Before the runs, the coefficients alpha, kappa, theta of the paper are
%   compared with offhooked_coeffs (our formulas): a check of the geometry.

    if nargin < 1, opts = struct(); end
    if ~isfield(opts, 'T'),     opts.T = 8; end
    if ~isfield(opts, 'cases'), opts.cases = 1:8; end
    setappdata(0, 'ctrl_phase', [0 0 0]);             % no phase shift
    addpath(fullfile(fileparts(mfilename('fullpath')), '..'));   % trailer2.m and params.m are one folder up
    Tr = @(name, varargin) trailer2(name, varargin{:});

    G(1).name = '(26)'; G(1).g = struct('d0', 0.1, 'd1', 0.1, 'l1', 1, 'l2', 1); G(1).c = [0.10495 0.0024 0.0112];
    G(2).name = '(27)'; G(2).g = struct('d0', 1,   'd1', 1,   'l1', 1, 'l2', 1); G(2).c = [1 1.5 4];
    z0_paper = [1; -2; 2; -1; 1];
    x0_small = [0.4; -0.4; 0.3; -0.2; 0.1];
    h_one = @(z) 1;
    h_log = @(z) abs(log(sum(z.^2)));

    for s = 1:2
        g = G(s).g;
        [a, k_, t_] = offhooked_coeffs(g.d0, g.d1, g.l1, g.l2);
        fprintf('geometry %s: alpha, kappa, theta  paper [%g %g %g]   our formulas [%.5g %.5g %.5g]\n', ...
                G(s).name, G(s).c, a, k_, t_);
    end

    ep = 0.1;  om = 2*pi/ep;  gam = 1;  k = [8 6 5];
    dt = min(0.02, 2*pi/(max(k)*om)/40);
    % case table: geometry index, start ('paper' | 'small'), h index (1 = one, 2 = log)
    C = [1 1 1; 1 1 2; 1 2 1; 1 2 2; 2 1 1; 2 1 2; 2 2 1; 2 2 2];
    startname = {'paper z0', 'small x0'};  hname = {'h = 1', 'h = |log|'};

    R = struct([]);
    for ci = opts.cases(:).'
        g = G(C(ci,1)).g;  c = G(C(ci,1)).c;
        if C(ci,2) == 1, x0 = Tr('z_to_x', z0_paper, g.d0, g.d1, g.l1, g.l2);
        else,            x0 = x0_small; end
        z0 = x_to_z(x0, g.d0, g.d1, g.l1, g.l2);  r0 = norm(z0);
        if C(ci,3) == 1, hf = h_one; else, hf = h_log; end
        label = sprintf('case %d: geometry %s, start %s, %s', ci, G(C(ci,1)).name, startname{C(ci,2)}, hname{C(ci,3)});
        fprintf('\n##### %s\n  x0 = [%s]   ||z0|| = %.4f\n', label, num2str(x0.', 4), r0);
        fprintf('  %6s %10s %10s\n', 't [s]', 'exact', 'approx');

        xc = x0(:);  zc = z0(:);  t_a = 0;  emin = 1;  umax = 0;  amax = 0;  verdict = 'not finished';
        while t_a < opts.T - 1e-9
            t_b = min(t_a + 1, opts.T);
            [~, xs, zs, u1] = Tr('simulate_case_physical', xc, g.d0, g.d1, g.l1, g.l2, c(1), c(2), c(3), ...
                                 gam, ep, om, k(1), k(2), k(3), t_a, t_b, dt, hf);
            [~, za]         = Tr('simulate_case', zc, c(1), c(2), c(3), gam, ep, om, k(1), k(2), k(3), ...
                                 t_a, t_b, dt, hf);
            re = sqrt(sum(zs.^2, 2))/r0;
            emin = min(emin, min(re));
            umax = max(umax, max(abs(u1)));
            amax = max(amax, max(max(abs(xs(:, 3:5)))));
            xc = xs(end, :).';  zc = za(end, :).';  t_a = t_b;
            fprintf('  %6.1f %10.4g %10.4g\n', t_a, re(end), norm(zc)/r0);
            if re(end) < 0.01, verdict = 'CONVERGED'; break, end
            if re(end) > 3,    verdict = 'DIVERGED';  break, end
        end
        fprintf('  -> %s: min exact %.4g, final exact %.4g, approx %.4g, max|u1| %.1f, largest |angle| %.0f deg\n', ...
                verdict, emin, re(end), norm(zc)/r0, umax, rad2deg(amax));
        R(end+1).case_id = ci;  R(end).label = label;  R(end).verdict = verdict;
        R(end).exact_min = emin;  R(end).exact_final = re(end);  R(end).approx_final = norm(zc)/r0;
    end
    fprintf('\nratios are ||z||/||z0||; "approx" = nilpotent approximation, "exact" = exact kinematics.\n');
end
