function P = params()
%PARAMS  Every number of the two-trailer study, in ONE place.
%
%   P = params()
%
%   Sections:  1 geometry   2 hardware limits   3 controller
%              4 Lyapunov function   5 scenario (start state, time grid)
%
%   Three rules that keep the rest of the code short and safe:
%     * Change the sampling period with   P = set_epsilon(P, eps)
%       and never by assigning P.epsilon alone -- P.omega = 2*pi/epsilon
%       would go stale (this has bitten before).
%     * Change the Lyapunov function by changing the lyapunov_h(...) line in
%       section 4. Nothing else has to move.
%     * P.x0 is in RADIANS: [x1; x2; phi; theta1; theta2].

    setup_paths();      % makes lib/ visible to everything that follows

    % =====================================================================
    %  1. GEOMETRY  [m]   (measured in the lab, 15 Sep 2026)
    % =====================================================================
    P.d0 = 0.115;
    P.d1 = 0.084;
    P.l1 = 0.082;
    P.l2 = 0.089;
    % for reference, the paper's simulation sets:
    %   set (26):  d0 = 0.1, d1 = 0.1, l1 = 1, l2 = 1
    %   set (27):  d0 = 1,   d1 = 1,   l1 = 1, l2 = 1

    % Coefficients of the nilpotent approximation (follow from geometry only)
    [P.alpha, P.kappa, P.theta, P.beta] = ...
        offhooked_coeffs(P.d0, P.d1, P.l1, P.l2);

    % =====================================================================
    %  2. HARDWARE LIMITS  (TurtleBot3 Burger)
    % =====================================================================
    P.u1_max = 0.22;       % driving velocity  [m/s]
    P.u2_max = 2.84;       % steering velocity [rad/s]

    % =====================================================================
    %  3. CONTROLLER
    % =====================================================================
    % Tuned for P.x0 below, V = gamma*(sqrt(r^2+d^2)-d), k = (14,5,1).
    % THREE things have to hold, and the third is the one that bites:
    %   (1) amplitude:   control_bound(z, P, gamma) <= (u1_max, u2_max)
    %   (2) contraction: eps*gamma*h(r) <= cap_product, i.e. gamma <=
    %                    cap*delta/eps  (h is bounded by 1/delta)
    %   (3) VALIDITY:    the nilpotent approximation has to still describe
    %                    the real kinematics at this eps.
    %
    % (3) is what caps eps. The controller dithers, and the excursion
    % inside one sampling interval grows with eps. Once it is no longer
    % small against the distance to the goal, the approximation stops
    % describing the plant -- and the plant is what the robot (and Gazebo)
    % actually is. Measured with the exact model (compare_physical.m),
    % 12 intervals:
    %   eps= 2 -> 0.03 m excursion, tracks to 3.3%
    %   eps= 3 -> 0.04 m           tracks to 5.4%   <- in use
    %   eps= 4 -> 0.06 m           tracks to 9.6%   (marginal, |u1| 0.222)
    %   eps= 6 -> 0.09 m           BREAKS
    %   eps=26 -> 0.43 m           diverges outright: ||z|| 0.84 -> 1.63,
    %             trailers wind up 23 turns, |u1| hits 0.358
    % So a large eps looks spectacular on the approximate model and is a
    % mirage. ALWAYS confirm a retune with compare_physical before
    % believing it.
    %
    % Full convergence to 0.1% of ||z0|| takes about 8.5 h of simulated
    % time at these settings. P.t_end below only shows the first stretch.
    P.gamma = 2.7498e-05;                 % gain of the fixed-gain controller
    P = set_epsilon(P, 3);                % sampling period [s]; sets P.omega too

    % frequency multipliers of the three oscillating terms
    % (was (18,4,1); (14,5,1) won the C1-admissible search at every
    %  epsilon tried -- see search_k.m)
    P.k12   = 14;
    P.k112  = 5;
    P.k1112 = 1;

    % Safety cap on eps*gamma*h for the ADAPTIVE-gain controller
    % (simulate(P, ..., 'adaptive'), compare_gain.m). The proof's linear
    % bound ||z(eps)|| <= |1-eps*gamma|*||z0|| predicts an exact hit at
    % eps*gamma = 1 and a breakdown around eps*gamma = 2.
    % 0.8 is because of convergence_map.m.
    P.cap_product = 0.8;

    % =====================================================================
    %  4. LYAPUNOV FUNCTION
    % =====================================================================
    % The controller only ever sees V through its gradient. For a radially
    % symmetric V = phi(r), r = ||z||, that collapses to one scalar factor:
    %
    %     a(z) = -G^-1(z)*grad V(z) = -gamma*h(z)*G^-1(z)*z
    %     h(z) = phi'(r) / (gamma*r)
    %
    % so CHANGING THE LYAPUNOV FUNCTION IS CHANGING h -- nothing else in
    % the code has to move. Options (see lyapunov_h.m):
    %
    %   lyapunov_h('quadratic')     V = gamma*r^2/2
    %   lyapunov_h('log', c)        V = gamma/(2c)*log(1 + c*r^2)   (Grushkovskaya's)
    %   lyapunov_h('sqrt', delta)   V = gamma*(sqrt(r^2 + delta^2) - delta)   <- now
    %
    % A NON-radial V (e.g. a weighted quadratic that fixes the very
    % different physical scales of z1..z5) does not fit this hook: there
    % grad V is no longer parallel to z and the closed forms in amps.m
    % must be re-derived with derive_symbolic.m.
    [P.h_fun, P.lyap_name] = lyapunov_h('sqrt', 0.005);

    % =====================================================================
    %  5. SCENARIO
    % =====================================================================
    P.t_start = 0;
    P.t_end   = 1800;     % NOT full convergence (that is ~8.5 h) -- just a
                          % viewable stretch; k12*omega = 29 rad/s here
    P.dt      = 2e-3;     % output grid / maximum solver step [s]

    % start state: [x1; x2; phi; theta1; theta2], RADIANS
    P.x0 = [-0.2; -0.2; deg2rad(45); deg2rad(45); deg2rad(45)];
    P.z0 = x_to_z(P.x0, P.d0, P.d1, P.l1, P.l2);
    % To start from a given z0 instead, overwrite P.z0 (and P.x0 =
    % z_to_x(P.z0, ...) if you also want the exact model to match).
end
