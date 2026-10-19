function P = params()
%PARAMS  Single source of truth: geometry, controller constants, scenario.
%   ONE-TRAILER (r = 3, 4-state) version. See params.m in "2 trailers/"
%   for the two-trailer (r = 4, 5-state) version -- do not have both
%   folders on the MATLAB path at once, the filenames collide.
%
%   P = params();
%   then unpack what you need, e.g.  alpha = P.alpha;
%
%   Rules that keep P consistent:
%     * change epsilon with  P = set_epsilon(P, eps)  -- it updates omega
%       too (P.omega is derived from P.epsilon below);
%     * change delta (or the Lyapunov function) in THIS file -- P.h_fun is
%       built from P.delta once, so editing P.delta afterwards does
%       nothing to h_fun;
%     * P.x0 is in DEGREES for the angles (x0_to_z0 takes degrees), unlike
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
    %
    % measured on the assembled hardware, 2026-09-15 (same car+hitch+trailer1
    % as in the "2 trailers" folder -- d1, l2 don't apply here, there is no
    % second trailer)
    P.d0 = 0.115;
    P.l1 = 0.082;

    P.alpha = one_trailer_coeff(P.d0, P.l1);

    % ---------------- controller --------------------------------------------
    % --- tuned for P.x0 below, V = gamma*(sqrt(r^2+d^2)-d), k=(3,2).
    % Largest gain meeting BOTH constraints at every point of the ray --
    % amplitude (control_bound <= 0.22, 2.84) and contraction
    % (eps*gamma*h(r) <= cap_product) -- and then VERIFIED by a full
    % oscillating simulation against the EXACT kinematics, run past the
    % target, not just to it:
    %   eps=3 -> gamma=1.39e-03, 0.1% of ||z0|| at 654 s
    %   eps=5 -> gamma=3.20e-03, 0.1% at 280 s, |u1| 0.214  <- in use
    %   eps=7 -> gamma=2.29e-03, 0.1% at 392 s
    % All of eps = 2..7 converge and stay inside the limits, so this is
    % not near a cliff. Re-check after ANY change to x0, k or delta.
    %   old values: gamma = 0.27, eps = 3.38
    P.gamma   = 3.2028e-03;
    P.epsilon = 5;
    P.omega   = 2*pi / P.epsilon;

    % Safety cap on eps*gamma*h for the ADAPTIVE-gain controller
    % (compare_gain.m, simulate_case_adaptive.m); see "2 trailers/params.m".
    P.cap_product = 0.8;

    % Only TWO frequency multipliers -- there is no third (fourth-order)
    % Lie bracket, so no k1112. The controller shape for r = 3 (the paper's
    % eq. (16) truncated after the second bracket) is confirmed. Which
    % non-resonance condition the pair (k12,k112) has to satisfy is NOT:
    % the paper's C1 is stated for three multipliers, and
    % heuristic_no_multiple.m is only a stand-in.
    % ---------------- Lyapunov function -------------------------------------
    % For a radially symmetric V = phi(r), r = ||z||, the whole function
    % enters the controller through one scalar factor:
    %     a(z) = -gamma*h(z)*G^-1(z)*z,   h(z) = phi'(r)/(gamma*r)
    % so CHANGING V IS CHANGING h_fun. Candidates:
    %   V = gamma*r^2/2              h = @(z) 1
    %   V = (gamma/2)*log(1+r^2)     h = @(z) 1/(1 + norm(z)^2)
    %   V = gamma*(sqrt(r^2+d^2)-d)  h = @(z) 1/sqrt(norm(z)^2 + d^2)   <- in use
    % The third one holds the control at the limit all the way in instead
    % of wasting headroom near the goal. For this start state, with the
    % largest admissible gain of each (gamma_ray.m) and the time from
    % eta_ray.m, to reach 0.1% of ||z0||: quadratic 27 min, log c=10 8 min,
    % sqrt delta=0.02 4.6 min. delta must be set BEFORE the handle is
    % created -- an anonymous function freezes the values it captures.
    P.delta = 0.02;
    P.h_fun = @(z) 1/sqrt(norm(z)^2 + P.delta^2);
    P.lyap_name = sprintf('sqrt, delta = %g', P.delta);   % label for printouts

    P.k12   = 3;      % search_k.m ranks the pairs that pass the heuristic
    P.k112  = 2;      % screen; (3,2) came out on top

    % ---------------- scenario ----------------------------------------------
    % Start state x = [x1; x2; phi; theta1] (goal = origin), angles in DEGREES
    % here. z0 is its image in privileged coordinates (x0_to_z0.m):
    %     z1 = x1,  z2 = phi,  z3 = -x2,
    %     z4 = l1*x2 - d0*l1^2/(d0+l1)*phi - l1^3/(d0+l1)*theta1
    %
    % Always sanity-check a new z0 by converting it back: an earlier default
    % [1;-2;2;-1] was just the 2-trailer z0 truncated to 4 entries and
    % corresponded to theta1 = 301.5 rad (48 full turns!).
    %
    % The gain above was tuned for THIS pose (0.2 m / 45 deg). Earlier
    % tunings in older notes belong to other poses, e.g. the gain 0.27 for
    % x = (0.036 m, 0.003 m, 2.1 deg, 1.4 deg); at the present pose it would
    % exceed |u1| by a factor of 5. Idea from that time, not tested since:
    % a start far along the track but with y, phi, theta1 near zero costs
    % little control, because a12 and a112 only contain z1 through products
    % with z2, z3.
    P.x0      = [-0.2; -0.2; 45; 45];   % [x1; x2; phi; theta1], angles in degrees
    P.z0      = x0_to_z0(P.x0, P.d0, P.l1);
    
    P.t_start = 0;
    P.t_end   = 400;   % 0.1% of ||z0|| is reached at ~280 s
    P.dt      = 2e-3;  % fastest oscillation is k12*omega = 3.8 rad/s now

    % ---------------- hardware limits ---------------------------------------
    P.u1_max = 0.22;      % TurtleBot3 Burger [m/s]
    P.u2_max = 2.84;      % TurtleBot3 Burger [rad/s]
end
