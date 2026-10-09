function z0 = x0_to_z0(x0, d0, d1, l1, l2)
%X0_TO_Z0  Convert a physical initial configuration into the privileged
%   (nilpotent-approximation) coordinates z0 used in params.m (P.z0).
%
%   x0 = [x1; x2; phi; theta1; theta2]   physical state
%        x1, x2   -- robot position [m], goal is the origin (0,0)
%        phi      -- robot heading [rad]
%        theta1   -- first trailer angle [rad]
%        theta2   -- second trailer angle [rad]
%
%   d0, d1, l1, l2  -- geometry [m], passed in directly (NOT read via
%                      params()) so this can safely be called FROM
%                      params.m itself without infinite recursion.
%
%   Example (inside params.m, after P.d0..P.l2 are already set):
%       P.x0 = [0.3; 0; deg2rad(45); deg2rad(45); deg2rad(45)];
%       P.z0 = x0_to_z0(P.x0, P.d0, P.d1, P.l1, P.l2);
%
%   Example (standalone, in the Command Window):
%       x0 = [0.3; 0; deg2rad(45); deg2rad(45); deg2rad(45)];
%       z0 = x0_to_z0(x0, 0.115, 0.084, 0.082, 0.089);

syms x1 x2 x3 x4 x5 real

f15 = sin(x3-x5)/l2 + (d1+l1)*sin(x4-x3)*cos(x4-x5)/(l1*l2);
f25 = d0*(d1+l1)*cos(x3-x4)*cos(x4-x5)/(l1*l2) - d0*cos(x3-x5)/l2;
f1  = [cos(x3); sin(x3); 0; sin(x3-x4)/l1; f15];
f2  = [0; 0; 1; -d0*cos(x3-x4)/l1; f25];

Xs  = [x1;x2;x3;x4;x5];
lie = @(f,g) simplify(jacobian(g,Xs)*f - jacobian(f,Xs)*g);

Y3 = lie(f1,f2);
Y4 = lie(f1,Y3);
Y5 = lie(f1,Y4);

M0 = double(subs([f1,f2,Y3,Y4,Y5], Xs, zeros(5,1)));

beta = d0*l2*(d0^2*(d1+l1)*(l1-l2-d1) + l1^3*(d0+l2)) / (6*(d0+l2)*l1^2);

x0 = x0(:);
y  = M0 \ x0;

z0    = y;
z0(5) = y(5) - beta*x0(3)^3;
end
