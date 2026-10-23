function [f1, f2] = fields_exact(x, d0, d1, l1, l2)
%FIELDS_EXACT  Exact kinematics of the car with two off-hooked trailers
%   (paper, eq. (2)-(3)):  x' = u1*f1(x) + u2*f2(x),
%   x = [x1; x2; phi; theta1; theta2], angles in radians.

    x3 = x(3);  x4 = x(4);  x5 = x(5);

    f15 = sin(x3 - x5)/l2 ...
          + (d1 + l1)*sin(x4 - x3)*cos(x4 - x5)/(l1*l2);
    f25 = d0*(d1 + l1)*cos(x3 - x4)*cos(x4 - x5)/(l1*l2) ...
          - d0*cos(x3 - x5)/l2;

    f1 = [cos(x3); sin(x3); 0; sin(x3 - x4)/l1; f15];
    f2 = [0; 0; 1; -d0*cos(x3 - x4)/l1; f25];
end
