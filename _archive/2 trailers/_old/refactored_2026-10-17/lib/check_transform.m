function check_transform(P)
%CHECK_TRANSFORM  Self-check of the coordinate transform x <-> z.
%
%   check_transform(P)  prints two errors and asserts both are tiny:
%     1. det(M0) must equal eq. (4) of the paper,
%            det F(0) = -(d0+l1)(d0+l2)(d1+l2) / (l1^4 l2^4),
%        i.e. the frame that defines the privileged coordinates is right.
%     2. z_to_x(x_to_z(x)) must give x back for a few test states.
%
%   x_to_z.m hard-codes the closed form of the frame, so this is the guard
%   against a typo there. Called by compare_physical.m and selftest.m.

    d0 = P.d0;  d1 = P.d1;  l1 = P.l1;  l2 = P.l2;

    [~, ~, ~, ~, M0] = offhooked_coeffs(d0, d1, l1, l2);
    det_paper = -(d0 + l1)*(d0 + l2)*(d1 + l2) / (l1^4*l2^4);
    err_det   = abs(det(M0) - det_paper) / abs(det_paper);

    % a few fixed, generic states (columns): position, heading, two angles
    X = [ 0.30 -0.50  0.10  0.70;
         -0.40  0.20  0.60 -0.30;
          0.50 -0.90  0.20  1.10;
         -0.70  0.40  0.90 -0.60;
          0.20  0.80 -1.00  0.50];
    err_rt = max(max(abs(z_to_x(x_to_z(X, d0, d1, l1, l2), d0, d1, l1, l2) - X)));

    fprintf('=== coordinate transform self-check ===\n');
    fprintf('  rel. error in det(M0) vs eq. (4)  : %.3e\n', err_det);
    fprintf('  max |z_to_x(x_to_z(x)) - x|       : %.3e\n', err_rt);
    assert(err_det < 1e-10, 'M0 does not match the determinant in the paper');
    assert(err_rt  < 1e-9,  'z_to_x is not the inverse of x_to_z');
    disp('   PASSED')
end
