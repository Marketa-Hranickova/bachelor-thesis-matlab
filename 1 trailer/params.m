function P = params(epsilon)
%PARAMS  Single source of truth: geometry, controller constants, scenario.
%   ONE-TRAILER (r = 3, 4-state) version. See params.m in "2 trailers/"
%   for the two-trailer (r = 4, 5-state) version -- do not have both
%   folders on the MATLAB path at once, the filenames collide.
%
%   P = params();         the numbers as set below
%   P = params(4);        the same, but with epsilon = 4 s (omega follows)
%   then unpack what you need, e.g.  alpha = P.alpha;
%
%   Rules that keep P consistent:
%     * change epsilon with  params(eps)  or by editing P.epsilon below --
%       P.omega is derived from P.epsilon, so never overwrite P.epsilon
%       on an existing P (P.omega would keep the old value and the
%       controller would oscillate at the wrong frequency, silently).
%       NOTE the gain is NOT retuned: P.gamma stays as tuned for the
%       epsilon below;
%     * change c (or delta, or the Lyapunov function) in THIS file -- P.h_fun is
%       built from P.c_log once, so editing P.c_log afterwards does
%       nothing to h_fun (and P.gamma has to be retuned for the new V);
%     * P.x0 is in DEGREES for the angles (x0_to_z0 below takes degrees), unlike
%       "2 trailers", where P.x0 is in radians. Convert with deg2rad before
%       calling the exact-kinematics simulator.
%
%   NOTE: alpha and gamma are BUILT-IN MATLAB functions. Returned as
%   struct fields on purpose -- see the note in the 2-trailer params.m
%   for why (bare names inside a function body resolve to the built-in
%   unless shadowed by an explicit assignment in that same function).

    % ---------------- geometry ---------------------------------------------
    %   d0 : car rear axle  -> hitch          [m]
    %   l1 : hitch          -> trailer axle   [m]
    % measured on the assembled hardware, 2026-09-15 (same car+hitch+trailer1
    % as in the "2 trailers" folder -- d1, l2 don't apply here, there is no
    % second trailer)
    P.d0 = 0.115;
    P.l1 = 0.082;

    P.alpha = one_trailer_coeff(P.d0, P.l1);   % local function at the end of this file

    % ---------------- controller --------------------------------------------
    % --- tuned for P.x0 below, V = gamma/(2c)*log(1+c*r^2), c = 10 (Victoria's), k=(3,2).
    % Largest gain meeting BOTH constraints at every point of the ray --
    % amplitude (control_bound <= 0.22, 2.84) and contraction
    % (eps*gamma*h(r) <= cap_product) -- and then VERIFIED by a full
    % oscillating simulation against the EXACT kinematics, run past the
    % target, not just to it. Re-check after ANY change to x0, k or the Lyapunov function.
    % Result for eps = 5:
    %   eps=5 -> gamma=2.1143e-02, 0.1% at 445 s (exact kinematics), |u1| 0.214, |u2| 0.37  <- in use
    %            (gradient-flow estimate 492 s; ray_gain gives gamma, eps*gamma*hmax = 0.11)
    P.gamma   = 2.1143e-02;
    P.epsilon = 5;
    if nargin >= 1 && ~isempty(epsilon), P.epsilon = epsilon; end   % params(eps)
    P.omega   = 2*pi / P.epsilon;

    % Safety cap on eps*gamma*h for the ADAPTIVE-gain controller
    % (compare_gain, simulate_case_adaptive); see "2 trailers/params.m".
    P.cap_product = 0.8;

    % Only TWO frequency multipliers -- there is no third (fourth-order)
    % Lie bracket, so no k1112. The controller shape for r = 3 (the paper's
    % eq. (16) truncated after the second bracket) is confirmed. Which
    % non-resonance condition the pair (k12,k112) has to satisfy is NOT:
    % the paper's C1 is stated for three multipliers, and
    % heuristic_no_multiple (in trailer1.m) is only a stand-in.
    % ---------------- Lyapunov function -------------------------------------
    % For a radially symmetric V = phi(r), r = ||z||, the whole function
    % enters the controller through one scalar factor:
    %     a(z) = -gamma*h(z)*G^-1(z)*z,   h(z) = phi'(r)/(gamma*r)
    % so CHANGING V IS CHANGING h_fun. Candidates:
    %   V = gamma*r^2/2              h = @(z) 1
    %   V = gamma/(2c)*log(1+c*r^2)  h = @(z) 1/(1 + c*norm(z)^2)   <- in use, c = 10 (Victoria's)
    %   V = gamma*(sqrt(r^2+d^2)-d)  h = @(z) 1/sqrt(norm(z)^2 + d^2)
    % Victoria's log V (her trailers2.m), c = 10. c must be set BEFORE the handle is created.
    P.c_log = 10;
    P.h_fun = @(z) 1/(1 + P.c_log*norm(z)^2);
    P.lyap_name = sprintf('log, c = %g (Victoria)', P.c_log);
    % the sqrt V used before (gamma = 3.2028e-03 for it):
    %   P.delta = 0.02;  P.h_fun = @(z) 1/sqrt(norm(z)^2 + P.delta^2);
    %   P.lyap_name = sprintf('sqrt, delta = %g', P.delta);

    P.k12   = 3;      % search_k ranks the pairs that pass the heuristic
    P.k112  = 2;      % screen; (3,2) came out on top

    % ---------------- scenario ----------------------------------------------
    % Start state x = [x1; x2; phi; theta1] (goal = origin), angles in DEGREES
    % here. z0 is its image in privileged coordinates (x0_to_z0, below):
    %     z1 = x1,  z2 = phi,  z3 = -x2,
    %     z4 = l1*x2 - d0*l1^2/(d0+l1)*phi - l1^3/(d0+l1)*theta1
    % Always sanity-check a new z0 by converting it back: an earlier default
    % [1;-2;2;-1] was just the 2-trailer z0 truncated to 4 entries and
    % corresponded to theta1 = 301.5 rad (48 full turns!).
    P.x0      = [-0.2; -0.2; 45; 45];   % [x1; x2; phi; theta1], angles in degrees
    P.z0      = x0_to_z0(P.x0, P.d0, P.l1);
    
    P.t_start = 0;
    P.t_end   = 400;
    P.dt      = 2e-3;

    % ---------------- where thesis_figures writes the PDFs ------------------
    P.fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');   % next to this file; created on demand

    % ---------------- hardware limits ---------------------------------------
    P.u1_max = 0.22;      % TurtleBot3 Burger [m/s]
    P.u2_max = 2.84;      % TurtleBot3 Burger [rad/s]
end

% ---- helpers (moved here from their own files; only params.m needs them) ----

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

function z0 = x0_to_z0(x0, d0, l1)
%X0_TO_Z0  Convert a physical initial configuration into the privileged
%   coordinates z0 used in params.m (P.z0). ONE-TRAILER (r=3, 4-state)
%   version -- see "2 trailers/x0_to_z0.m" for the two-trailer one.
%
%   x0 = [x1; x2; phi_deg; theta1_deg]   physical state
%        x1, x2              -- robot position [m], goal is the origin (0,0)
%        phi_deg             -- robot heading, IN DEGREES
%        theta1_deg          -- trailer angle, IN DEGREES
%
%   d0, l1  -- geometry [m], passed in directly (NOT read via params())
%              so this can safely be called FROM params.m itself without
%              infinite recursion.
%
%   Closed form (verified independently, matches the comment already in
%   params.m and the from-scratch Lie-bracket derivation -- angles
%   converted to radians internally before use):
%       z1 = x1
%       z2 = phi
%       z3 = -x2
%       z4 = l1*x2 - d0*l1^2/(d0+l1)*phi - l1^3/(d0+l1)*theta1
%   No nonlinear correction is needed here (unlike the 2-trailer z5) --
%   the degree of nonholonomy is only r = 3, see DERIVATION_4STATE.md.
%
%   Example (inside params.m, after P.d0, P.l1 are already set):
%       P.x0 = [0.15; 0; 0; 0];       % 0 deg, 0 deg
%       P.z0 = x0_to_z0(P.x0, P.d0, P.l1);
%
%   Example (standalone, in the Command Window):
%       x0 = [0.2; 0.2; 45; 45];      % 45 deg, 45 deg
%       z0 = x0_to_z0(x0, 0.115, 0.082);

x0 = x0(:);
x1 = x0(1); x2 = x0(2);
x3 = deg2rad(x0(3));
x4 = deg2rad(x0(4));

z0 = [ x1;
       x3;
      -x2;
       l1*x2 - d0*l1^2/(d0+l1)*x3 - l1^3/(d0+l1)*x4 ];
end
