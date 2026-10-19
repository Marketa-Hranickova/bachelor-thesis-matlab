function [f1, f2] = fields_exact(x, d0, l1)
%FIELDS_EXACT  Exact kinematics of the car with ONE off-hooked trailer:
%   x' = u1*f1(x) + u2*f2(x),  x = [x1; x2; phi; theta1], angles in radians.
%   These are rows 1..4 of the 5-state two-trailer model (the first four
%   rows never involve the second trailer; see derive_report.m).

    x3 = x(3);  x4 = x(4);

    f1 = [cos(x3); sin(x3); 0;  sin(x3 - x4)/l1];
    f2 = [0;       0;       1; -d0*cos(x3 - x4)/l1];
end
