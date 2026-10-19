function varargout = derive(task, varargin)
%DERIVE  Symbolic derivation for the ONE-TRAILER system. Needs the Symbolic
%   Math Toolbox; none of the simulations need it (they are in trailer1.m).
%
%   derive('symbolic')   Lie brackets, controllability, amplitudes a(z) of the
%                        4-state system, checked against amps in trailer1.m
%   derive('report')     derive_privileged on BOTH systems (4 and 5 state),
%                        checked against the paper and params.m (a few minutes)
%   R = derive('privileged', f1, f2, X)
%                        the Bellaiche construction itself, for any driftless
%                        two-input system (see the help of derive_privileged below)
%
%   Nothing here was changed in content when the files were merged; the three
%   used to be derive_symbolic.m, derive_report.m and derive_privileged.m.

    if nargin < 1 || isempty(task)
        help derive
        return
    end
    fns = struct('symbolic', @derive_symbolic, 'report', @derive_report, ...
                 'privileged', @derive_privileged);
    if ~isfield(fns, task)
        error('derive:unknownTask', 'Unknown task "%s" -- use symbolic, report or privileged.', task);
    end
    f = fns.(task);
    if nargout == 0
        f(varargin{:});
    else
        [varargout{1:nargout}] = f(varargin{:});
    end
end


function derive_symbolic()
% DERIVE_SYMBOLIC  Lie brackets, controllability and the amplitudes a(z) of
% the off-hooked ONE-TRAILER system (r = 3, 4 states) -- symbolically.
%
%   Needs the Symbolic Math Toolbox and takes a while (simplify). Run it to
%   see or re-check the algebra behind the controller; simulations do not
%   need it. (This was PART 1 and 2 of calculations.m, before everything moved into
%   trailer1.m.)
%
%   Derivation: DERIVATION_4STATE.md in the thesis notes (same method as
%   the paper -- Bellaiche/Jean/Risler privileged coordinates -- redone for
%   this 4-state system, since the paper itself only treats r = 4).
%
%   Part 1  vector fields g1, g2 of the nilpotent approximation, their Lie
%           brackets, the matrix G(z) = [g1, g2, [g1,g2], [g1,[g1,g2]]] and
%           det G (= 1).
%   Part 2  a(z) = -G^-1(z)*grad Q(z) for the quadratic Q = gamma/2*||z||^2,
%           checked against amps in trailer1.m. Any radial V gives the same formulas
%           times the scalar h(z) (see params.m).

    P = params();
    alpha = P.alpha;  gamma = P.gamma;

    fprintf('geometry : d0=%.3f  l1=%.3f  [m]\n', P.d0, P.l1);
    fprintf('coeff    : alpha=%.6f\n', alpha);

    % =========================================================================
    %  PART 1 -- symbolic system, Lie brackets, controllability
    % =========================================================================
    syms z1 z2 z3 z4 real
    syms alpha_s real
    syms gamma_s positive

    z_sym = [z1; z2; z3; z4];

    % Vector fields of the nilpotent approximation (r = 3). The system is
    %     z' = u1*g1(z) + u2*g2(z)
    % One state fewer than the 2-trailer case -- no z5, no kappa, no theta.
    % (numerically the same fields as fields_approx in trailer1.m)
    g1_sym = [1; 0; -z2; -z3];
    g2_sym = [0; 1; 0; -alpha_s*z2^2];

    fprintf('\n=== Vector field g1(z) ===\n'); disp(g1_sym)
    fprintf('\n=== Vector field g2(z) ===\n'); disp(g2_sym)

    dg1dz = jacobian(g1_sym, z_sym);
    dg2dz = jacobian(g2_sym, z_sym);

    % [g1, g2]
    Lie_g1g2_sym = simplify(dg2dz*g1_sym - dg1dz*g2_sym);
    fprintf('\n=== [g1, g2] ===\n'); disp(Lie_g1g2_sym)

    % [g1, [g1, g2]] -- this is the LAST bracket needed for r = 3. There is no
    % [g1,[g1,[g1,g2]]] here (that belongs to a 5th state that doesn't exist).
    dLie_g1g2_dz    = jacobian(Lie_g1g2_sym, z_sym);
    Lie_g1_g1g2_sym = simplify(dLie_g1g2_dz*g1_sym - dg1dz*Lie_g1g2_sym);
    fprintf('\n=== [g1, [g1, g2]] ===\n'); disp(Lie_g1_g1g2_sym)

    % The four brackets span the tangent space -> the system is controllable.
    G_sym = [g1_sym, g2_sym, Lie_g1g2_sym, Lie_g1_g1g2_sym];
    fprintf('\n=== G matrix ===\n'); disp(G_sym)

    detG = simplify(det(G_sym));
    fprintf('\n=== det(G) ===\n'); disp(detG)
    fprintf('(DERIVATION_4STATE.md says this should be exactly 1, for all z)\n');

    % =========================================================================
    %  PART 2 -- the amplitudes a(z), for h = 1
    % =========================================================================
    % Lyapunov function Q(z) = gamma/2 * ||z||^2, now over 4 states.
    Q_sym     = gamma_s/2 * (z1^2 + z2^2 + z3^2 + z4^2);
    gradQ_sym = jacobian(Q_sym, z_sym).';
    fprintf('\n=== grad Q(z) ===\n'); disp(gradQ_sym)

    G_inv_sym = simplify(inv(G_sym));

    % a(z) = -h * G^{-1} * grad Q,  with h = 1
    a_sym = simplify(-G_inv_sym*gradQ_sym);

    fprintf('\n=== a(z) components, h = 1 ===\n')
    labels = {'a1   ', 'a2   ', 'a12  ', 'a112 '};
    for i = 1:4
        fprintf('%s = ', labels{i}); disp(a_sym(i))
    end

    % ------------------------------------------------------------------ check
    % Every simulation uses the closed forms in amps (trailer1.m). Check them against the
    % symbolic result at a generic point, so a typo in either place is caught.
    z_test = [0.37; -0.21; 0.44; -0.13];

    a_from_sym = double(subs(a_sym, ...
        [z_sym.', gamma_s, alpha_s], ...
        [z_test.', gamma, alpha]));

    % amps with h = 1, to match the quadratic V derived above
    [b1, b2, b12, b112] = trailer1('amps', z_test, alpha, gamma, @(z) 1);
    a_from_amps = [b1; b2; b12; b112];

    res = max(abs(a_from_sym - a_from_amps));
    fprintf('\n=== symbolic vs amps: max |difference| = %.3e ===\n', res);
    assert(res < 1e-10, 'a(z) mismatch -- symbolic result and amps disagree');
    disp('   PASSED')

    fprintf('\n=== Done. ===\n')
end


function derive_report()
% DERIVE_REPORT
% Runs derive_privileged on BOTH systems and checks the output against
% what is independently known.
%
% The engine is told only f1 and f2. It is never told r, the weights, the
% correction, alpha or beta -- it derives them. This script then compares:
%
%   5-state : against the paper (arXiv:2403.03341), whose beta and alpha
%             are published. Reproducing beta from scratch is the real test.
%   4-state : against the alpha in params.m and amps in trailer1.m, i.e.
%             against the code the rest of this folder runs on.
%
% Runtime: a few minutes, most of it the 5-state case.

    syms x1 x2 x3 x4 x5 real

    % measured geometry
    d0 = sym(115)/1000;  d1 = sym(84)/1000;
    l1 = sym(82)/1000;   l2 = sym(89)/1000;

    f15 = sin(x3-x5)/l2 + (d1+l1)*sin(x4-x3)*cos(x4-x5)/(l1*l2);
    f25 = d0*(d1+l1)*cos(x3-x4)*cos(x4-x5)/(l1*l2) - d0*cos(x3-x5)/l2;
    F1  = [cos(x3); sin(x3); 0; sin(x3-x4)/l1; f15];
    F2  = [0;       0;       1; -d0*cos(x3-x4)/l1; f25];

    %% ========================= 4-STATE =====================================
    fprintf('================= 4-STATE (one trailer) =================\n');
    R4 = derive_privileged(F1(1:4), F2(1:4), [x1 x2 x3 x4]);
    assert(R4.ok, 'weighted-order test failed');
    assert(isequal(R4.r, 3) && isequal(R4.growth, [2 3 4]), 'r / growth wrong');
    assert(isequal(simplify(R4.corr), sym(zeros(4,1))), ...
           'a correction came out non-zero for the 4-state case');

    fprintf('\n  derived g1 = '); disp(R4.g1.');
    fprintf('  derived g2 = '); disp(R4.g2.');

    alpha4 = simplify(-R4.g2(4)/R4.z(2)^2);
    fprintf('  derived alpha = '); disp(alpha4);
    assert(isAlways(simplify(alpha4 - d0*(d0+l1)/2) == 0), 'alpha wrong');
    fprintf('  == d0*(d0+l1)/2                              -> OK\n');
    Pp = params();
    assert(abs(double(alpha4) - Pp.alpha) < 1e-12);
    fprintf('  alpha in params.m agrees                     -> OK\n');

    % amplitudes against amps (trailer1.m)
    zt = [0.31; -0.17; 0.42; -0.08];
    ref = double(subs(R4.a, [R4.z.' R4.gam], [zt.' 1]));
    try
        [b1,b2,b12,b112] = trailer1('amps', zt, double(alpha4), 1, @(z) 1);
        assert(max(abs([b1;b2;b12;b112] - ref)) < 1e-12, 'amps disagrees');
        fprintf('  amps agrees numerically                      -> OK\n');
    catch ME
        fprintf(2,'  amps check FAILED: %s\n', ME.message);
        fprintf(2,'  (if the message is about params or trailer1: the 2-trailer\n');
        fprintf(2,'   folder is shadowing this one -- run "which -all params")\n');
    end

    %% ========================= 5-STATE =====================================
    fprintf('\n================= 5-STATE (two trailers) =================\n');
    fprintf('(slower -- this is the one that must reproduce the paper''s beta)\n');
    R5 = derive_privileged(F1, F2, [x1 x2 x3 x4 x5]);
    assert(R5.ok, 'weighted-order test failed');
    assert(isequal(R5.r, 4) && isequal(R5.growth, [2 3 4 5]), 'r / growth wrong');

    % only z5 may carry a correction
    assert(isequal(simplify(R5.corr(1:4)), sym(zeros(4,1))), ...
           'an unexpected correction appeared on z1..z4');
    beta_derived = simplify(R5.corr(5)/R5.z(2)^3);
    beta_paper   = d0*l2*(d0^2*(d1+l1)*(l1-l2-d1) + l1^3*(d0+l2))/(6*(d0+l2)*l1^2);

    fprintf('\n  derived beta = '); disp(beta_derived);
    fprintf('  paper''s beta = '); disp(simplify(beta_paper));
    assert(isAlways(simplify(beta_derived - beta_paper) == 0), ...
           'DERIVED BETA DOES NOT MATCH THE PAPER');
    fprintf('  exact match                                  -> OK\n');
    fprintf('  (z5 = y5 - beta*z2^3, and z2 = x3, so this is eq.(9) of the paper)\n');

    alpha5 = simplify(-R5.g2(4)/R5.z(2)^2);
    fprintf('\n  derived alpha (5-state) = '); disp(alpha5);
    fprintf('  numerically %.10g\n', double(alpha5));

    fprintf('\n==========================================================\n');
    fprintf('  BOTH SYSTEMS DERIVED FROM f1, f2 ALONE.\n');
    fprintf('  The 5-state run reproduced the published beta.\n');
    fprintf('==========================================================\n');
end


function R = derive_privileged(f1, f2, X, verbose)
%DERIVE_PRIVILEGED  Bellaiche construction for a driftless two-input system.
%
%   R = derive_privileged(f1, f2, X)
%
%   Given the exact vector fields of  x' = u1*f1(x) + u2*f2(x)  this DERIVES,
%   without being told the answer:
%
%     R.r        degree of nonholonomy
%     R.growth   growth vector
%     R.w        weights of the privileged coordinates
%     R.F0       the frame [f1, f2, [f1,f2], ...] at the origin
%     R.detF0    its determinant
%     R.z        the privileged coordinate symbols z1..zn
%     R.corr     the NONLINEAR correction of each coordinate,
%                z_q = y_q - R.corr(q)   (zero where none is needed)
%     R.g1,R.g2  vector fields of the nilpotent approximation
%     R.a        amplitudes a(z) = -gamma*G^-1*z  (R.gam is the symbol)
%     R.ok       true if every coordinate passed the weighted-order test
%
%   Nothing above is hard-coded: run it on a different model and you get
%   that model's answer. Reproduces the paper's beta for the 5-state car
%   with two off-hooked trailers (see derive('report')).
%
%   Method: Bellaiche's recursion (Bellaiche/Jean/Risler, Geometry of
%   Nonholonomic Systems). The fields are Taylor-truncated in x before the
%   weighted grading, which is what keeps the 5-state case tractable.

    if nargin < 4, verbose = true; end
    n  = numel(X);  Xc = X(:);
    O  = num2cell(zeros(1,n));
    at0 = @(M) subs(M, X(:).', zeros(1,n));
    lie = @(f,g) expand(jacobian(g,Xc)*f - jacobian(f,Xc)*g);

    % ---- frame: bracket with f1 until the rank at the origin reaches n ---
    Y = {f1(:), f2(:)};  len = [1 1];  B = f2(:);
    while rank(at0([Y{:}])) < n
        B = lie(f1(:), B);  Y{end+1} = B;  len(end+1) = len(end)+1; %#ok<AGROW>
        assert(numel(Y) <= n+2, 'frame never spans -- is the system controllable?');
    end
    r = max(len);  w = len;
    growth = arrayfun(@(s) rank(at0([Y{len<=s}])), 1:r);

    F0    = at0([Y{:}]);
    detF0 = simplify(det(F0));
    assert(~isAlways(detF0 == 0), 'det F(0) = 0 -- not controllable at the origin');
    Finv  = simplify(inv(F0));

    if verbose
        fprintf('  n = %d,  r = %d,  growth = (%s),  weights = (%s)\n', ...
                n, r, join(string(growth),','), join(string(w),','));
        fprintf('  det F(0) = '); disp(detF0);
    end

    % ---- polynomialise in x (total degree <= r+1) before any grading -----
    P1 = sym(zeros(n,1));  P2 = sym(zeros(n,1));
    for i = 1:n
        P1(i) = poly_trunc(f1(i), Xc, r+1);
        P2(i) = poly_trunc(f2(i), Xc, r+1);
    end
    Yp = cell(size(Y));
    for j = 1:numel(Y)
        V = sym(zeros(n,1));
        for i = 1:n, V(i) = poly_trunc(Y{j}(i), Xc, r+1); end
        Yp{j} = V;
    end

    % ---- step 1: linear coordinates, step 2: Bellaiche correction --------
    y = expand(Finv*Xc);
    z = sym('z', [n 1], 'real');
    Ld = @(f,V) expand(jacobian(f, Xc)*V);

    corr = sym(zeros(n,1));
    for q = 1:n
        A = multi_indices(w(1:q-1), w(q));
        acc = sym(0);
        for i = 1:size(A,1)
            a = A(i,:);
            g = y(q);
            for j = q-1:-1:1                 % rightmost operator acts first
                for k = 1:a(j), g = Ld(g, Yp{j}); end
            end
            c = simplify(subs(g, X(:).', zeros(1,n)));
            if isAlways(c == 0, 'Unknown', 'false'), continue; end
            acc = acc + c/prod(factorial(a))*prod(z(1:q-1).'.^a);
        end
        corr(q) = expand(acc);
        if verbose && ~isequal(corr(q), sym(0))
            fprintf('  CORRECTION  z%d = y%d - (', q, q); 
            fprintf('%s)\n', char(simplify(corr(q))));
        end
    end
    if verbose && isequal(simplify(corr), sym(zeros(n,1)))
        fprintf('  no nonlinear correction needed\n');
    end

    % ---- zdot in z -------------------------------------------------------
    x_of_z = expand(F0*(z + corr));
    syms u1 u2 real
    Jz  = eye(n) - jacobian(corr, z);
    raw = expand(Finv*(P1*u1 + P2*u2));
    zdot = expand(Jz*subs(raw, Xc.', x_of_z.'));

    % ---- leading weighted term of each component -------------------------
    syms ep positive
    lead = sym(zeros(n,1));  ok = true;
    if verbose, fprintf('  weighted-order test:\n'); end
    for q = 1:n
        e = expand(subs(zdot(q), z.', (ep.^w(:)).'.*z.'));
        low = NaN;
        for j = 0:(r+2)
            if j == 0, cj = simplify(subs(e, ep, 0));
            else,      cj = simplify(subs(diff(e, ep, j), ep, 0)/factorial(j)); end
            if ~isAlways(cj == 0, 'Unknown', 'false')
                low = j;  lead(q) = cj;  break
            end
        end
        pass = (low == w(q)-1);  ok = ok && pass;
        if verbose
            if pass, tag = 'OK'; else, tag = '*** FAIL ***'; end
            fprintf('    zdot%d: lowest %d, required %d  %s\n', ...
                    q, low, w(q)-1, tag);
        end
    end

    % ---- nilpotent approximation and amplitudes --------------------------
    g1 = simplify(subs(lead, [u1 u2], [1 0]));
    g2 = simplify(subs(lead, [u1 u2], [0 1]));
    Gz = sym(zeros(n));  Gz(:,1) = g1;  Gz(:,2) = g2;  Bz = g2;
    for j = 3:n
        Bz = expand(jacobian(Bz,z)*g1 - jacobian(g1,z)*Bz);
        Gz(:,j) = simplify(Bz);
    end
    syms gam positive
    a = simplify(-gam*(Gz\z));

    R = struct('r',r,'growth',growth,'w',w,'F0',F0,'detF0',detF0, ...
               'z',z,'corr',corr,'g1',g1,'g2',g2,'G',Gz,'a',a, ...
               'gam',gam,'ok',ok,'lead',lead);
end


function e = poly_trunc(e, Xc, deg)
%POLY_TRUNC  multivariate Taylor of e about x = 0, total degree <= deg
    syms t_ positive
    e = subs(e, Xc.', (t_*Xc).');
    e = expand(subs(taylor(e, t_, 0, 'Order', deg+1), t_, 1));
end


function A = multi_indices(wlow, wmax)
%MULTI_INDICES  all non-zero integer tuples a with  sum(a.*wlow) < wmax
    q = numel(wlow);
    A = zeros(0, q);
    if q == 0, return; end
    caps = floor((wmax-1)./wlow);
    idx  = zeros(1, q);
    while true
        if any(idx) && sum(idx.*wlow) < wmax, A(end+1,:) = idx; end %#ok<AGROW>
        k = q;
        while k >= 1
            idx(k) = idx(k) + 1;
            if idx(k) <= caps(k), break; end
            idx(k) = 0;  k = k - 1;
        end
        if k < 1, break; end
    end
end
