function P = set_epsilon(P, epsilon)
%SET_EPSILON  Change the sampling period AND the base frequency together.
%
%   P = set_epsilon(P, 4)
%
%   omega = 2*pi/epsilon is derived from epsilon. Writing P.epsilon = ...
%   directly leaves P.omega at the old value, and the controller then
%   oscillates at the wrong frequency without any error. Use this instead.

    P.epsilon = epsilon;
    P.omega   = 2*pi / epsilon;
end
