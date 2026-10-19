function [alpha, kappa, theta, beta, M0] = offhooked_coeffs(d0, d1, l1, l2)
%OFFHOOKED_COEFFS  Everything that follows from the geometry alone.
%
%   [alpha, kappa, theta, beta, M0] = offhooked_coeffs(d0, d1, l1, l2)
%
%   alpha, kappa, theta : coefficients of the nilpotent approximation
%   beta                : cubic correction in z5
%   M0                  : the frame F(0) = [f1, f2, [f1,f2], [f1,[f1,f2]],
%                         [f1,[f1,[f1,f2]]]] at the origin, which defines
%                         the privileged coordinates (see x_to_z.m, z_to_x.m)
%
%   Asking for only four outputs (as trailers2.m does) works too.
%   The closed forms are from Zuyev & Grushkovskaya (arXiv:2403.03341).

    alpha = ( d0^3*( l1^3 + l1^2*l2 - l1*l2^2 - d1*l2*(d1 + l2) ) ...
            + d0*l1^3*( d0*(l1 + 2*l2) + l2*(l1 + l2) ) ) ...
            / ( 2*(d0 + l2)*l1^3 );

    kappa = d0^2*( 1.5*l2^2*(d0 + l1)*(d1 + l1) ...
            - l2*( l1^3 + 0.5*d0*l1^2 - d1*l1*(d0 + 1.5*d1) ...
                   - 1.5*d0*d1^2 ) ...
            + 0.5*d0*d1^2*l1 ) ...
            / ( 3*(d0 + l2)*l1^3 );

    theta = d0^2*( l2^2*(d1 + l1)*(d0 + l1) ...
            + ( d0*d1 + (d0 + d1)*l1 )*d1*l2 ...
            + d0*d1^2*l1 ) ...
            / ( (d0 + l2)*l1^3 );

    beta  = ( d0*l2*( d0^2*(d1 + l1)*(l1 - l2 - d1) ...
            + l1^3*(d0 + l2) ) ) ...
            / ( 6*(d0 + l2)*l1^2 );

    % ---- frame at the origin ---------------------------------------------
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
