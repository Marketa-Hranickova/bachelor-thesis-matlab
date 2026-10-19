function [f1, f2] = fields_exact(x, P)
%FIELDS_EXACT  Vector fields of the EXACT kinematics at physical state x.
%
%   [f1, f2] = fields_exact(x, P)
%
%   The real plant (car with two off-hooked trailers, paper eq. (2)-(3)):
%       x' = u1*f1(x) + u2*f2(x),    x = [x1; x2; phi; theta1; theta2]
%   with angles in RADIANS. This is what the robot (and Gazebo) really
%   does; fields_approx.m is the truncated model the controller is
%   designed on. The gap between the two is what compare_physical.m shows.

    x3 = x(3);  x4 = x(4);  x5 = x(5);

    f15 = sin(x3 - x5)/P.l2 ...
          + (P.d1 + P.l1)*sin(x4 - x3)*cos(x4 - x5)/(P.l1*P.l2);
    f25 = P.d0*(P.d1 + P.l1)*cos(x3 - x4)*cos(x4 - x5)/(P.l1*P.l2) ...
          - P.d0*cos(x3 - x5)/P.l2;

    f1 = [cos(x3); sin(x3); 0; sin(x3 - x4)/P.l1; f15];
    f2 = [0; 0; 1; -P.d0*cos(x3 - x4)/P.l1; f25];
end
