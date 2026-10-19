function x = z_to_x(z, d0, l1)
%Z_TO_X  Privileged coordinates -> physical state. Exact inverse of
%   x_to_z.m (angles in RADIANS). Useful for plotting where the vehicle
%   actually drives, in metres and radians, instead of the abstract z.

    z = z(:);
    A = d0*l1^2/(d0+l1);
    B = l1^3/(d0+l1);
    x = [ z(1);
         -z(3);
          z(2);
          (-l1*z(3) - A*z(2) - z(4))/B ];
end
