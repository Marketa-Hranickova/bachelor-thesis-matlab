function P = params()
%PARAMS  Single source of truth: geometry, controller constants, scenario.
%   ONE-TRAILER (r = 3, 4-state) version. See params.m in "2 trailers/"
%   for the two-trailer (r = 4, 5-state) version -- do not have both
%   folders on the MATLAB path at once, the filenames collide.
%
%   P = params();
%   then unpack what you need, e.g.  alpha = P.alpha;
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

    % Safety cap on eps*gamma, for an adaptive-gain controller later (see
    % "2 trailers/compare_gain.m" for what this is for).
    P.cap_product = 0.8;

    % Only TWO frequency multipliers now -- there is no third (fourth-order)
    % Lie bracket, so no k1112. Reused the 2-trailer folder's (k12,k112)
    % pair as a starting point; whether the article's condition C1 even
    % applies in this 2-multiplier form is UNCONFIRMED (see calculations.m).
    % ---------------- Lyapunov function -------------------------------------
    % For a radially symmetric V = phi(r), r = ||z||, the whole function
    % enters the controller through one scalar factor:
    %     a(z) = -gamma*h(z)*G^-1(z)*z,   h(z) = phi'(r)/(gamma*r)
    % so CHANGING V IS CHANGING h_fun. Candidates:
    %   V = gamma*r^2/2              h = @(z) 1
    %   V = (gamma/2)*log(1+r^2)     h = @(z) 1/(1 + norm(z)^2)
    %   V = gamma*(sqrt(r^2+d^2)-d)  h = @(z) 1/sqrt(norm(z)^2 + d^2)   <- in use
    % The third one holds the control at the limit all the way in instead
    % of wasting headroom near the goal; measured ~4.6x faster than the
    % plain quadratic. delta must be set BEFORE the handle is created --
    % an anonymous function freezes the values it captures.
    P.delta = 0.02;
    P.h_fun = @(z) 1/sqrt(norm(z)^2 + P.delta^2);

    P.k12   = 3;      % a C1-free search over admissible pairs picked (3,2)
    P.k112  = 2;      % independently -- see the heuristic in search_k.m

    % ---------------- scenario ----------------------------------------------
    % z0 corresponds to a PHYSICAL starting pose x = (xC0,yC0,phi,theta1) =
    % (0.036 m, 0.003 m, 0.036 rad = 2.1 deg, 0.025 rad = 1.4 deg), converted
    % forward through the privileged-coordinate map z1=x1, z2=x3, z3=-x2,
    % z4 = l1*x2 - d0*l1^2/(d0+l1)*x3 - l1^3/(d0+l1)*x4.
    % IMPORTANT: the earlier default [1;-2;2;-1] was just the 2-trailer
    % z0 truncated to 4 entries, never checked against physical reality --
    % it corresponds to theta1 = 301.5 rad (48 full turns!). Always sanity
    % check a new z0 this way before trusting "how far from the goal" talk.
    %
    % UPDATE: pushed further along-track while staying in limits. z1
    % (along-track distance) barely costs anything IF heading and lateral
    % offset are both near zero -- a12 has a z1*z2 cross term and a112 has
    % a z1*z3 cross term, so z1 alone only feeds the *linear* a1 term. So
    % instead of shrinking x uniformly, x = (x1, y, phi, theta1) =
    % (0.15 m, 0.5 mm, 0, 0) keeps y/phi/theta1 essentially perfect (0.5 mm
    % lateral offset is inside the ~0.05 mm mocap noise floor territory,
    % see RESEARCH_LOG.md) and lets x1 be 4x farther than the old 3.6 cm
    % point. Verified by full simulation (not just the initial-state
    % algebraic bound, which is NOT strictly conservative across the whole
    % run -- see finding from 17 Sep).
    %
    % STALE-NOTE WARNING: the verification that used to be recorded here
    % (max|u1| = 0.186, inside the limit) belongs to the OLD pose
    % (0.036 m, 0.003 m, 2.1 deg, 1.4 deg) described above, NOT to the
    % P.x0 now in force. At the current 0.2 m / 45 deg pose the old gain
    % 0.27 gives a bound of |u1| = 1.12, i.e. 5.1x over. The gain above
    % was retuned for THIS pose.
    P.x0      = [-0.2; -0.2; 45; 45];   % [x1; x2; phi; theta1], goal = origin
    P.z0      = x0_to_z0(P.x0, P.d0, P.l1);
    
    P.t_start = 0;
    P.t_end   = 400;   % 0.1% of ||z0|| is reached at ~280 s
    P.dt      = 2e-3;  % fastest oscillation is k12*omega = 3.8 rad/s now

    % ---------------- hardware limits ---------------------------------------
    P.u1_max = 0.22;      % TurtleBot3 Burger [m/s]
    P.u2_max = 2.84;      % TurtleBot3 Burger [rad/s]
end
