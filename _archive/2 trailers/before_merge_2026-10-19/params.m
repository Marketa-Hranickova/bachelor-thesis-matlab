function P = params()
%PARAMS  Every number of the two-trailer study, in one struct P.
%
%   P = params()   is called by every script and simulator; change numbers
%   HERE and nowhere else.
%
%   Two rules that keep P consistent:
%     * change epsilon with  P = set_epsilon(P, eps)  -- it updates omega
%       too (P.omega is derived from P.epsilon below);
%     * change delta (or the Lyapunov function) in THIS file -- P.h_fun is
%       built from P.delta once, so editing P.delta afterwards does
%       nothing to h_fun.

    % ---------------- geometry [m], measured 15 Sep 2026 --------------------
    P.d0 = 0.115;
    P.d1 = 0.084;
    P.l1 = 0.082;
    P.l2 = 0.089;
    % paper set (26):  0.1, 0.1, 1, 1
    % paper set (27):  1,   1,   1, 1

    [P.alpha, P.kappa, P.theta, P.beta] = ...
        offhooked_coeffs(P.d0, P.d1, P.l1, P.l2);

    % ---------------- controller --------------------------------------------
    % --- tuned for P.x0 below, V = gamma*(sqrt(r^2+d^2)-d), k=(14,5,1).
    % THREE things have to hold, and the third is the one that bites:
    %   (1) amplitude:   control_bound(z, gamma*h) <= (u1_max, u2_max)
    %   (2) contraction: eps*gamma*h(r) <= cap_product, i.e. gamma <=
    %                    cap*delta/eps  (h is bounded by 1/delta)
    %   (3) VALIDITY:    the nilpotent approximation has to still describe
    %                    the real kinematics at this eps.
    %
    % (3) is what caps eps. The controller dithers, and the excursion
    % inside one sampling interval grows with eps. Once it is no longer
    % small against the distance to the goal, the approximation stops
    % describing the plant -- and the plant is what the robot (and Gazebo)
    % actually is. Measured against simulate_case_physical, 12 intervals:
    %   eps= 2 -> 0.03 m excursion, tracks to 3.3%
    %   eps= 3 -> 0.04 m           tracks to 5.4%   <- in use
    %   eps= 4 -> 0.06 m           tracks to 9.6%   (marginal, |u1| 0.222)
    %   eps= 6 -> 0.09 m           BREAKS
    %   eps=26 -> 0.43 m           diverges outright: ||z|| 0.84 -> 1.63,
    %             trailers wind up 23 turns, |u1| hits 0.358
    % So a large eps looks spectacular on simulate_case and is a mirage.
    % ALWAYS confirm a retune with compare_physical before believing it.
    %
    % Full convergence to 0.1% of ||z0|| takes about 8.5 h of simulated
    % time at these settings. P.t_end below only shows the first stretch.
    P.gamma   = 0.5;
    P.epsilon = 2;
    P.omega   = 2*pi / P.epsilon;

    % Safety cap on eps*gamma*h for the ADAPTIVE-gain controller
    % (compare_gain.m, simulate_case_adaptive.m) and compare_V.m. The proof's linear
    % bound ||z(eps)|| <= |1-eps*gamma|*||z0|| predicts an exact hit at
    % eps*gamma = 1 and a breakdown around eps*gamma = 2.
    % 0.8 = safety margin below 1, read off convergence_map.m.
    P.cap_product = 0.8;

    % ---------------- Lyapunov function -------------------------------------
    % The controller only ever sees V through its gradient. For a radially
    % symmetric V = phi(r) with r = ||z|| that collapses to one scalar
    % factor:
    %
    %     a(z) = -G^-1(z)*grad V(z) = -gamma*h(z)*G^-1(z)*z
    %     h(z) = phi'(r) / (gamma*r)
    %
    % so CHANGING THE LYAPUNOV FUNCTION IS CHANGING h_fun -- nothing else
    % in the code has to move. Candidates:
    %
    %   V = gamma*r^2/2                 h = @(z) 1                 (quadratic)
    %   V = gamma/(2c)*log(1+c*r^2)     h = @(z) 1/(1 + c*norm(z)^2)  (log, Grushkovskaya)
    %   V = gamma*(sqrt(r^2+d^2)-d)     h = @(z) 1/sqrt(norm(z)^2 + d^2) (sqrt)  <- in use
    %
    % compare_V.m ranks these under the velocity limits.
    %
    % A NON-radial V (e.g. a weighted quadratic that fixes the very
    % different physical scales of z1..z5) does not fit this hook: there
    % grad V is no longer parallel to z and the closed forms in amps.m
    % must be re-derived (derive_symbolic.m).

    %P.h_fun = @(z) 1;
    P.delta = 0.005;  % with eps=3, (2) is slack, so delta can be small
    P.h_fun = @(z) 1/sqrt(norm(z)^2 + P.delta^2);
    P.lyap_name = sprintf('sqrt, delta = %g', P.delta);   % label for printouts

    P.k12   = 14;      % was (18,4,1); (14,5,1) won the C1-admissible
    P.k112  = 5;       % search at every epsilon tried
    P.k1112 = 1;

    % ---------------- scenario ----------------------------------------------
    P.t_start = 0;
    P.t_end   = 1800;  % NOT full convergence (that is ~8.5 h) -- just a
                       % viewable stretch; k12*omega = 29 rad/s here
    P.dt      = 2e-3;
    % start state [x1 x2 phi theta1 theta2], angles in RADIANS
    P.x0 = [-0.2; -0.2; deg2rad(45); deg2rad(45); deg2rad(45)];
    P.z0 = x_to_z(P.x0, P.d0, P.d1, P.l1, P.l2);

    % ---------------- hardware limits, TurtleBot3 Burger --------------------
    P.u1_max = 0.22;      % m/s
    P.u2_max = 2.84;      % rad/s
end
