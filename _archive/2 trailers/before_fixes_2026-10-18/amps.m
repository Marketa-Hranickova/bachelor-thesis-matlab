function [a1, a2, a12, a112, a1112] = amps(zf, alpha, kappa, theta, gamma_val, h_fun)
    h     = h_fun(zf);
    a1    = -gamma_val*h*zf(1);
    a2    = -gamma_val*h*zf(2);
    a12   = -gamma_val*h*(zf(1)*zf(2) + zf(3));
    a112  = -gamma_val*h*(zf(1)*zf(3) + zf(4) + alpha*zf(2)^3);
    a1112 = -gamma_val*h*( (alpha + 2*kappa - theta)*zf(1)*zf(2)^3 ...
                         + (alpha + 3*kappa)*zf(2)^2*zf(3) ...
                         + zf(1)*zf(4) + zf(5) );
end
