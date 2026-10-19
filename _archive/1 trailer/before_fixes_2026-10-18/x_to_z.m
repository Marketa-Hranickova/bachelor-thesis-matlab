function z = x_to_z(x, d0, l1)
%X_TO_Z  Physical state -> privileged coordinates. ONE-TRAILER (r=3).
%
%   x = [x1; x2; phi; theta1]   angles in RADIANS (unlike x0_to_z0.m,
%                               which takes degrees, for convenience when
%                               setting a scenario up by hand)
%
%   This is the version to call inside a loop -- the controller speaks z,
%   the robot measures x, so this runs at the start of every sampling
%   interval. The controls u1, u2 need NO conversion: they are the
%   physical driving / steering velocities and are identical in both
%   descriptions, because x -> z changes only the STATE coordinates.

    x = x(:);
    z = [ x(1);
          x(3);
         -x(2);
          l1*x(2) - d0*l1^2/(d0+l1)*x(3) - l1^3/(d0+l1)*x(4) ];
end
