function M0 = privileged_frame(d0, d1, l1, l2)
%PRIVILEGED_FRAME  The matrix F(0) = [f1, f2, [f1,f2], [f1,[f1,f2]],
%   [f1,[f1,[f1,f2]]]] of the exact system, evaluated at the origin.
%   This is the frame that defines the privileged coordinates:
%       y = F(0) \ x        (paper, eq. (5)-(6))
%
%   Closed form, derived once from the symbolic Lie brackets and then
%   hard-coded, so that no Symbolic Toolbox call is needed at run time.
%   Cross-checked against x0_to_z0.m in compare_physical.m.
%
%   Sanity check on the algebra:
%       det(M0) = -(d0+l1)*(d0+l2)*(d1+l2) / (l1^4*l2^4)
%   which is exactly eq. (4) of the paper.

    M0 = zeros(5);

    M0(1,1) =  1;
    M0(2,3) = -1;
    M0(3,2) =  1;

    M0(4,2) = -d0/l1;
    M0(4,3) = -(d0 + l1)/l1^2;
    M0(4,4) = -(d0 + l1)/l1^3;
    M0(4,5) = -(d0 + l1)/l1^4;

    % Row 5 follows the pattern  [ d0*d1*S_k + l1*l2*(d0+d1)*S_{k-1}
    %                              + l1^2*l2^2*S_{k-2} ] / (l1^j*l2^j)
    % with S_k = l1^k + l1^(k-1)*l2 + ... + l2^k.
    M0(5,2) = d0*d1 / (l1*l2);
    M0(5,3) = ( d0*d1*(l1 + l2) + l1*l2*(d0 + d1) ) / (l1^2*l2^2);
    M0(5,4) = ( d0*d1*(l1^2 + l1*l2 + l2^2) ...
                + l1*l2*(d0 + d1)*(l1 + l2) ...
                + l1^2*l2^2 ) / (l1^3*l2^3);
    M0(5,5) = ( d0*d1*(l1^3 + l1^2*l2 + l1*l2^2 + l2^3) ...
                + l1*l2*(d0 + d1)*(l1^2 + l1*l2 + l2^2) ...
                + l1^2*l2^2*(l1 + l2) ) / (l1^4*l2^4);
end
