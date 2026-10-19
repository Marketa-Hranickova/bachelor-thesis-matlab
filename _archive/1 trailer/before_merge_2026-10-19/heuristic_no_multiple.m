function ok = heuristic_no_multiple(a, b)
%HEURISTIC_NO_MULTIPLE  Practical placeholder non-resonance check for the
%   2-multiplier (r=3, one-trailer) reduction: neither k12 nor k112 may be
%   an integer multiple of the other.
%
%   THIS IS NOT THE PAPER'S CONDITION C1. That condition is proven only
%   for the 5-state (r=4), 3-multiplier construction -- see the caveat in
%   params.m in this folder. Which condition the 2-multiplier case needs
%   is OPEN (to ask Grushkovskaya); this is just a common-sense stand-in.
%   Do not present results filtered by this as more than a heuristic screen.

    ok = mod(a, b) ~= 0 && mod(b, a) ~= 0;
end
