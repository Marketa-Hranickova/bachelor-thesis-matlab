function [u1_val, u2_val] = u_ctrl(t, zf, omega, k12, k112, k1112, ...
                                   alpha, kappa, theta, gamma_val, h_fun)
%U_CTRL  Both control components at time t, for the state zf frozen at the
%   start of the current sampling interval.
%
%   u1 is the driving velocity [m/s], u2 the steering velocity [rad/s].
%   They are PHYSICAL inputs: the same u1, u2 drive the approximate
%   z-system and the exact x-system. Only the STATE needs converting.
%
%   Both components share the same five amplitudes, so they are computed
%   together -- that also halves the work inside the ode45 right-hand side
%   compared with the old split u1_ctrl / u2_ctrl pair.

    [a1, a2, a12, a112, a1112] = amps(zf, alpha, kappa, theta, gamma_val, h_fun);

    T12   = (2*k12*omega*abs(a12))^(1/2);
    T112  = (2*k112*omega)^(2/3) * abs(a112)^(1/3);
    T1112 = (k1112*omega)^(3/4) * abs(a1112)^(1/4);

    u1_val = a1 ...
        + T12   * cos(k12*omega*t) ...
        + T112  * cos(k112*omega*t) ...
        + 2*T1112 * cos(k1112*omega*t);

    u2_val = a2 ...
        + T12   * sign(a12)   * sin(k12*omega*t) ...
        - 2*T112  * sign(a112)  * cos(2*k112*omega*t) ...
        - 6*T1112 * sign(a1112) * sin(3*k1112*omega*t);
end
