function [u1, u2] = u_ctrl(t, c, P)
%U_CTRL  Control inputs at time t, for the amplitudes frozen in c.
%
%   [u1, u2] = u_ctrl(t, c, P)
%
%   c = ctrl_frame(z_j, P, gamma) holds the numbers frozen at the start of
%   the current sampling interval (see ctrl_frame.m); t is absolute time.
%
%   u1 is the driving velocity [m/s], u2 the steering velocity [rad/s].
%   They are PHYSICAL inputs: the same u1, u2 drive the approximate
%   z-system and the exact x-system. Only the STATE needs converting.

    w = P.omega;

    u1 = c.a1 ...
        + c.T12   * cos(P.k12*w*t) ...
        + c.T112  * cos(P.k112*w*t) ...
        + 2*c.T1112 * cos(P.k1112*w*t);

    u2 = c.a2 ...
        + c.T12   * c.s12   * sin(P.k12*w*t) ...
        - 2*c.T112  * c.s112  * cos(2*P.k112*w*t) ...
        - 6*c.T1112 * c.s1112 * sin(3*P.k1112*w*t);
end
