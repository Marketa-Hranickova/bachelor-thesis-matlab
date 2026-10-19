function ok = heuristic_no_multiple(a, b)
%HEURISTIC_NO_MULTIPLE  Practical placeholder non-resonance check for the
%   2-multiplier (r=3, one-trailer) reduction: neither k12 nor k112 may be
%   an integer multiple of the other.
%
%   THIS IS NOT THE PAPER'S CONDITION C1. That condition is proven only
%   for the 5-state (r=4), 3-multiplier construction -- see the caveat in
%   params.m and calculations.m in this folder. Whether ANY condition (and
%   if so, which one) is required for this 2-multiplier case is
%   UNCONFIRMED; this is just a common-sense stand-in until Grushkovska
%   confirms the real one. Do not present results filtered by this as
%   more than a preliminary/heuristic screen.

    ok = mod(a, b) ~= 0 && mod(b, a) ~= 0;
end
