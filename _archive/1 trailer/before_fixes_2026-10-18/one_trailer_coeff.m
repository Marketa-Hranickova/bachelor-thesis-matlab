function alpha = one_trailer_coeff(d0, l1)
%ONE_TRAILER_COEFF  The single coupling coefficient for the CAR + ONE
%   off-hooked trailer system (r = 3, 4 privileged coordinates).
%
%   alpha = d0*(d0+l1)/2
%
%   Only d0 (car rear axle -> hitch) and l1 (hitch -> trailer axle) are
%   needed -- there is no second hitch/trailer, so no d1, l2.
%
%   This is a DIFFERENT, SIMPLER formula than offhooked_coeffs.m in the
%   "2 trailers" folder -- it is not that formula with d1=l2=0, it comes
%   from redoing the privileged-coordinates derivation for the 4-state
%   system from scratch (see AIOS/.../DERIVATION_4STATE.md).
%
%   Sanity checks (see DERIVATION_4STATE.md):
%     alpha -> 0 as d0 -> 0, and the system then reduces exactly to
%     chained form (z4' = -z3*u1) -- an independent confirmation that the
%     derivation is right.
%     alpha > 0 always (no sign-changing degenerate geometry like the
%     2-trailer alpha has at d0 ~ 227 mm).

    alpha = d0*(d0 + l1)/2;
end
