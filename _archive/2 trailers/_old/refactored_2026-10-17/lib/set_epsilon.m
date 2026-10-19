function P = set_epsilon(P, epsilon)
%SET_EPSILON  Change the sampling period AND the base frequency together.
%
%   P = set_epsilon(P, 2)    ->  P.epsilon = 2,  P.omega = 2*pi/2
%
%   Use this instead of assigning P.epsilon: omega must always equal
%   2*pi/epsilon, and forgetting to update it is a silent error.

    P.epsilon = epsilon;
    P.omega   = 2*pi/epsilon;
end
