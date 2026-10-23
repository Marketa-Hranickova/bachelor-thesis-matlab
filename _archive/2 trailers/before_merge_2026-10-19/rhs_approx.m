function dz = rhs_approx(t, z, zf, g1_fun, g2_fun, omega, k12, k112, k1112, ...
                         alpha, kappa, theta, gamma_val, h_fun)
%RHS_APPROX  Closed-loop right-hand side of the nilpotent approximation,
%   z' = u1*g1(z) + u2*g2(z), with u1, u2 from u_ctrl at the frozen state zf.
%   Used by simulate_case.m and simulate_case_adaptive.m.
%
%   g1_fun, g2_fun may be [] -- then the fields come from fields_approx.m
%   (the normal case). Passing your own handles still works, for
%   experiments with modified fields.

    [u1, u2] = u_ctrl(t, zf, omega, k12, k112, k1112, ...
                      alpha, kappa, theta, gamma_val, h_fun);
    if isempty(g1_fun)
        [g1, g2] = fields_approx(z, alpha, kappa, theta);
    else
        g1 = g1_fun(z, alpha, kappa, theta);
        g2 = g2_fun(z, alpha, kappa, theta);
    end
    dz = g1*u1 + g2*u2;
end
