function z = x_to_z(x, d0, d1, l1, l2)
%X_TO_Z  Physical state  ->  privileged (nilpotent) coordinates.
%
%   z = x_to_z(x, d0, d1, l1, l2)
%
%   x = [x1; x2; phi; theta1; theta2]   physical state, angles in RADIANS
%   z = [z1; z2; z3; z4; z5]            privileged coordinates
%
%   x may also be a 5-by-N matrix (one state per column); z then is 5-by-N.
%
%   Two steps (Bellaiche's construction):
%       y    = M0 \ x                   linear part, M0 = frame F(0)
%       z5   = y5 - beta*phi^3          the one nonlinear correction
%   This is pure arithmetic, so it is fast enough to run at the start of
%   every sampling interval -- THIS is the version to port to the robot:
%   the controller needs z, the robot measures x.
%
%   The controls u1, u2 themselves need NO conversion -- they are the
%   physical driving / steering velocities and are identical in both
%   descriptions, because x -> z is a change of STATE coordinates only.
%
%   (The slow symbolic derivation this was validated against lives in
%   _derivation/x0_to_z0.m.)

    [~, ~, ~, beta, M0] = offhooked_coeffs(d0, d1, l1, l2);

    if isvector(x), x = x(:); end

    z        = M0 \ x;
    z(5, :)  = z(5, :) - beta*x(3, :).^3;
end
