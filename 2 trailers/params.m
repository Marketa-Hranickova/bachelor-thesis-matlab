function P = params()


    P.d0 = 0.115;
    P.d1 = 0.084;
    P.l1 = 0.082;
    P.l2 = 0.089;
    % paper set (26):  0.1, 0.1, 1, 1
    % paper set (27):  1,   1,   1, 1

    [P.alpha, P.kappa, P.theta, P.beta] = ...
        offhooked_coeffs(P.d0, P.d1, P.l1, P.l2);

    % ---------------- controller --------------------------------------------
    P.gamma   = 0.45;
    P.epsilon = 2;
    P.omega   = 2*pi / P.epsilon;

    % Safety cap on eps*gamma for the ADAPTIVE-gain controller
    % (compare_gain.m, simulate_case_adaptive.m). The proof's linear
    % bound ||z(eps)|| <= |1-eps*gamma|*||z0|| predicts an exact hit at
    % eps*gamma = 1 and a breakdown around eps*gamma = 2.
    % 0.8 is because of convergence map
    P.cap_product = 0.8;

    P.k12   = 18;
    P.k112  = 4;
    P.k1112 = 1;

    % ---------------- scenario ----------------------------------------------
    %P.z0      = [0.2;0.7854;0.2;0.01697;-0.00043];
    P.t_start = 0;
    P.t_end   = 30;
    P.dt      = 1e-4;
    P.x0 = [-0.2; -0.2; deg2rad(45); deg2rad(45); deg2rad(45)];
    P.z0 = x_to_z(P.x0, P.d0, P.d1, P.l1, P.l2);

    % ---------------- hardware limits for urtleBot3 burger ---------------------------------------
    P.u1_max = 0.22;      % m/s
    P.u2_max = 2.84;      % rad/s
end
