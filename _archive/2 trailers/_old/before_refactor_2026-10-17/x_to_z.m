function z = x_to_z(x, d0, d1, l1, l2)
%X_TO_Z  Physical state  ->  privileged (nilpotent) coordinates.
%
%   x = [x1; x2; phi; theta1; theta2]   physical state, angles in RADIANS
%   z = [z1; z2; z3; z4; z5]            privileged coordinates
%
%   Same map as x0_to_z0.m, but in closed form: x0_to_z0.m re-derives the
%   Lie brackets symbolically on every call (fine once, far too slow in a
%   loop), while this one is pure arithmetic. THIS is the version to port
%   to the robot: the controller needs z, the robot measures x, so this
%   runs at the start of every sampling interval.
%
%   The controls u1, u2 themselves need NO conversion -- they are the
%   physical driving / steering velocities and are identical in both
%   descriptions, because x -> z is a change of STATE coordinates only.

    [~, ~, ~, beta, M0] = offhooked_coeffs(d0, d1, l1, l2);

    x = x(:);
    y = M0 \ x;

    z    = y;
    z(5) = y(5) - beta*x(3)^3;
end
