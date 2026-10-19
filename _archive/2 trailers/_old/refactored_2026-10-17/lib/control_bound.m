function [b1, b2, terms1, terms2] = control_bound(z, P, gamma)
%CONTROL_BOUND  Worst-case |u1|, |u2| at state z (triangle inequality).
%
%   [b1, b2, terms1, terms2] = control_bound(z, P, gamma)
%
%   The bound is the constant part plus the magnitude of every oscillating
%   term -- all cosines/sines at their peak at once, so it is conservative,
%   but the real peaks come close to it (see plots / compare_physical).
%
%   terms1 = [ |a1|, T12,   T112,   2*T1112 ]    b1 = sum(terms1)
%   terms2 = [ |a2|, T12, 2*T112, 6*T1112 ]      b2 = sum(terms2)
%
%   Uses P.h_fun, so the bound is correct for whatever Lyapunov function is
%   set in params.m. It scales with the same frozen numbers the controller
%   uses (ctrl_frame.m), so the two cannot drift apart.

    c = ctrl_frame(z, P, gamma);

    terms1 = [abs(c.a1), c.T12,   c.T112, 2*c.T1112];
    terms2 = [abs(c.a2), c.T12, 2*c.T112, 6*c.T1112];

    b1 = sum(terms1);
    b2 = sum(terms2);
end
