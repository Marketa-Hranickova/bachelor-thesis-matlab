function [a1, a2, a12, a112, a1112] = controller_amps(z, alpha, kappa, theta, gamma_val)


    a1    = -gamma_val*z(1);
    a2    = -gamma_val*z(2);
    a12   = -gamma_val*(z(1)*z(2) + z(3));
    a112  = -gamma_val*(z(1)*z(3) + z(4) + alpha*z(2)^3);
    a1112 = -gamma_val*( (alpha + 2*kappa - theta)*z(1)*z(2)^3 ...
                       + (alpha + 3*kappa)*z(2)^2*z(3) ...
                       + z(1)*z(4) + z(5) );
end
