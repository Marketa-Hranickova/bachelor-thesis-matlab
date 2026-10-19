function x = z_to_x(z, d0, d1, l1, l2)
%Z_TO_X  Privileged coordinates  ->  physical state. Exact inverse of x_to_z.
%
%   x = z_to_x(z, d0, d1, l1, l2)
%
%   Useful for plotting where the vehicle actually drives (metres and
%   radians) instead of the abstract z coordinates. z may be 5-by-N.
%
%   The inverse is explicit because row 3 of M0 is [0 1 0 0 0], so the
%   heading is simply  phi = x(3) = z(2)  and the cubic correction can be
%   undone without solving anything.

    [~, ~, ~, beta, M0] = offhooked_coeffs(d0, d1, l1, l2);

    if isvector(z), z = z(:); end

    y        = z;
    y(5, :)  = z(5, :) + beta*z(2, :).^3;      % phi = z(2), see above

    x = M0*y;
end
