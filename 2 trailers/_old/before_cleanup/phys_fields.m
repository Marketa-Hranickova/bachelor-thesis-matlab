function [f1, f2] = phys_fields(x, d0, d1, l1, l2)
%PHYS_FIELDS  The two vector fields of the EXACT kinematics of the car
%   with two off-hooked trailers (paper, eq. (2)-(3)):
%
%       x' = u1*f1(x) + u2*f2(x),   x = [x1; x2; phi; theta1; theta2]
%
%   u1 is the driving velocity, u2 the steering velocity of the front
%   axle. These are the same fields x0_to_z0.m builds symbolically; here
%   they are plain numbers so they can be fed to ode45.
%
%   Note the contrast with g1/g2 in calculations.m: those are the
%   NILPOTENT APPROXIMATION of these fields near the origin. Everything
%   the controller is designed on lives there; this file is the real
%   plant.

    x3 = x(3);  x4 = x(4);  x5 = x(5);

    f15 = sin(x3 - x5)/l2 ...
          + (d1 + l1)*sin(x4 - x3)*cos(x4 - x5)/(l1*l2);
    f25 = d0*(d1 + l1)*cos(x3 - x4)*cos(x4 - x5)/(l1*l2) ...
          - d0*cos(x3 - x5)/l2;

    f1 = [cos(x3); sin(x3); 0; sin(x3 - x4)/l1; f15];
    f2 = [0; 0; 1; -d0*cos(x3 - x4)/l1; f25];
end
