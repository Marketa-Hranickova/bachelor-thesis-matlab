function [a1, a2, a12, a112, a1112] = amps(z, P, gamma)
%AMPS  The five controller amplitudes a(z) = -h(z) * G^-1(z) * grad Q(z).
%
%   [a1, a2, a12, a112, a1112] = amps(z, P, gamma)
%
%   z      state in privileged coordinates (5x1)
%   P      parameter struct (uses P.alpha, P.kappa, P.theta, P.h_fun)
%   gamma  gain
%
%   What the subscripts mean: each amplitude belongs to one direction the
%   robot can move in, and the subscript spells out the Lie bracket that
%   produces it.
%       a1, a2   -> driving / steering directly            f1, f2
%       a12      -> sideways "parallel parking" motion     [f1, f2]
%       a112     ->                                        [f1,[f1,f2]]
%       a1112    ->                                        [f1,[f1,[f1,f2]]]
%   a = -G^-1 grad Q says: "to descend the Lyapunov function, here is how
%   much of each available direction to use".
%
%   The closed forms below are the symbolic result for h = 1; the
%   state-dependent factor h(z) scales all five equally (it enters only
%   through the gradient of a radial V). derive_symbolic.m re-derives
%   them from the vector fields and checks they agree with this file.

    h     = P.h_fun(z);
    g     = gamma*h;                                   % effective gain here

    a1    = -g*z(1);
    a2    = -g*z(2);
    a12   = -g*(z(1)*z(2) + z(3));
    a112  = -g*(z(1)*z(3) + z(4) + P.alpha*z(2)^3);
    a1112 = -g*( (P.alpha + 2*P.kappa - P.theta)*z(1)*z(2)^3 ...
               + (P.alpha + 3*P.kappa)*z(2)^2*z(3) ...
               + z(1)*z(4) + z(5) );
end
