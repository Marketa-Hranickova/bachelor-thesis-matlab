function z0 = x0_to_z0(x0, d0, l1)
%X0_TO_Z0  Convert a physical initial configuration into the privileged
%   coordinates z0 used in params.m (P.z0). ONE-TRAILER (r=3, 4-state)
%   version -- see "2 trailers/x0_to_z0.m" for the two-trailer one.
%
%   x0 = [x1; x2; phi_deg; theta1_deg]   physical state
%        x1, x2              -- robot position [m], goal is the origin (0,0)
%        phi_deg             -- robot heading, IN DEGREES
%        theta1_deg          -- trailer angle, IN DEGREES
%
%   d0, l1  -- geometry [m], passed in directly (NOT read via params())
%              so this can safely be called FROM params.m itself without
%              infinite recursion.
%
%   Closed form (verified independently, matches the comment already in
%   params.m and the from-scratch Lie-bracket derivation -- angles
%   converted to radians internally before use):
%       z1 = x1
%       z2 = phi
%       z3 = -x2
%       z4 = l1*x2 - d0*l1^2/(d0+l1)*phi - l1^3/(d0+l1)*theta1
%   No nonlinear correction is needed here (unlike the 2-trailer z5) --
%   the degree of nonholonomy is only r = 3, see DERIVATION_4STATE.md.
%
%   Example (inside params.m, after P.d0, P.l1 are already set):
%       P.x0 = [0.15; 0; 0; 0];       % 0 deg, 0 deg
%       P.z0 = x0_to_z0(P.x0, P.d0, P.l1);
%
%   Example (standalone, in the Command Window):
%       x0 = [0.2; 0.2; 45; 45];      % 45 deg, 45 deg
%       z0 = x0_to_z0(x0, 0.115, 0.082);

x0 = x0(:);
x1 = x0(1); x2 = x0(2);
x3 = deg2rad(x0(3));   % phi
x4 = deg2rad(x0(4));   % theta1

z0 = [ x1;
       x3;
      -x2;
       l1*x2 - d0*l1^2/(d0+l1)*x3 - l1^3/(d0+l1)*x4 ];
end
