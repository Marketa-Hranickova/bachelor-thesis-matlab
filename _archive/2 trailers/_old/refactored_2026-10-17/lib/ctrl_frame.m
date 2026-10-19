function c = ctrl_frame(zf, P, gamma)
%CTRL_FRAME  Everything the controller freezes at the start of an interval.
%
%   c = ctrl_frame(zf, P, gamma)
%
%   The controller (Zuyev & Grushkovskaya, Theorem 1) works in sampling
%   intervals of length epsilon. At the start of each interval it looks at
%   the state zf ONCE, computes the five amplitudes a(zf), and then holds
%   them while the oscillating terms run:
%
%       u1 = a1 + T12*cos(k12*w*t) + T112*cos(k112*w*t) + 2*T1112*cos(k1112*w*t)
%       u2 = a2 + T12*s12*sin(k12*w*t) - 2*T112*s112*cos(2*k112*w*t)
%                                      - 6*T1112*s1112*sin(3*k1112*w*t)
%
%   This function returns those frozen numbers in a struct c:
%       c.a1, c.a2                  constant parts
%       c.a12, c.a112, c.a1112      bracket amplitudes
%       c.T12, c.T112, c.T1112      sizes of the oscillating terms
%       c.s12, c.s112, c.s1112      signs of a12, a112, a1112
%
%   The size of a term for a bracket of depth d is
%       T ~ (const * k * w)^((d-1)/d) * |a|^(1/d)
%   so deep brackets cost much more frequency (w^(3/4) for the deepest).
%
%   u_ctrl(t, c, P) evaluates the controls inside the interval;
%   control_bound(z, P, gamma) turns c into a worst-case bound.

    [c.a1, c.a2, c.a12, c.a112, c.a1112] = amps(zf, P, gamma);

    c.T12   = (2*P.k12*P.omega*abs(c.a12))^(1/2);
    c.T112  = (2*P.k112*P.omega)^(2/3) * abs(c.a112)^(1/3);
    c.T1112 = (P.k1112*P.omega)^(3/4)  * abs(c.a1112)^(1/4);

    c.s12   = sign(c.a12);
    c.s112  = sign(c.a112);
    c.s1112 = sign(c.a1112);
end
