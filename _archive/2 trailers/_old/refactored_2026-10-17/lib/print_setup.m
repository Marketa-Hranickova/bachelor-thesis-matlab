function print_setup(P)
%PRINT_SETUP  One-screen summary of the current parameters.
%
%   print_setup(P)  -- every run script starts with this, so each printout
%   (and each pasted output in an e-mail) says exactly what was simulated.

    fprintf('geometry  : d0=%.3f  d1=%.3f  l1=%.3f  l2=%.3f  [m]\n', ...
            P.d0, P.d1, P.l1, P.l2);
    fprintf('coeffs    : alpha=%.6f  kappa=%.6f  theta=%.6f\n', ...
            P.alpha, P.kappa, P.theta);
    fprintf('controller: epsilon=%.2f s (omega=%.4f rad/s)  gamma=%.4g  k=(%d,%d,%d)\n', ...
            P.epsilon, P.omega, P.gamma, P.k12, P.k112, P.k1112);
    fprintf('Lyapunov  : V = %s\n', P.lyap_name);
    fprintf('start     : x0=[%.3f %.3f m  %.1f %.1f %.1f deg],  ||z0||=%.4f\n', ...
            P.x0(1), P.x0(2), rad2deg(P.x0(3)), rad2deg(P.x0(4)), rad2deg(P.x0(5)), ...
            norm(P.z0));
    fprintf('horizon   : t in [%g, %g] s,  dt = %g s\n', P.t_start, P.t_end, P.dt);
    fprintf('limits    : |u1| <= %.2f m/s,  |u2| <= %.2f rad/s\n', P.u1_max, P.u2_max);
end
