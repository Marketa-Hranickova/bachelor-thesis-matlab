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
%   with two off-hooked trailers (see derive_report.m).
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
