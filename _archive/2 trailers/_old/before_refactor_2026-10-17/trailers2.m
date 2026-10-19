%% ================================================================
% two_trailer_loglyap.m
%
% Nilpotent approximation + oscillatory Lie-bracket controller
%
% Lyapunov function:
%
%   Q_c(z) = gamma/(2*c) * log(1 + c*||z||^2)
%
% Gradient:
%
%   grad Q = gamma*z / (1 + c*||z||^2)
%
% Hence the paper amplitudes are multiplied by
%
%   s(z) = 1/(1 + c*||z||^2)
%
% Larger c:
%   -> smaller control amplitudes
%   -> slower convergence
%
% Required external helper functions:
%
%   offhooked_coeffs.m
%   x0_to_z0.m
%
% ================================================================

clear;
clc;
close all;


%% ================================================================
% 1. PHYSICAL PARAMETERS
% ================================================================

d0 = 0.115;
d1 = 0.084;
l1 = 0.082;
l2 = 0.089;


[alpha,kappa,theta,beta_geom] = ...
    offhooked_coeffs(d0,d1,l1,l2);


%% ================================================================
% 2. CONTROLLER PARAMETERS
% ================================================================

gamma = 0.1;

epsilon = 2.0;

omega = 2*pi/epsilon;


k12   = 8;
k112  = 6;
k1112 = 5;


%% ================================================================
% 3. LYAPUNOV PARAMETER
%
% Q_c(z) = gamma/(2*c) log(1+c||z||^2)
%
% Start with c = 10.
%
% If controls are too large:
%
%   20
%   50
%   100
%   200
%   ...
%
% ================================================================

c_lyap = 10;


%% ================================================================
% 4. INITIAL CONDITION
% ================================================================

x0 = [ ...
    0.2;
    0.2;
    deg2rad(45);
    deg2rad(45);
    deg2rad(45) ];


z0 = x_to_z( ...
    x0,d0,d1,l1,l2);

z0 = z0(:);


%% ================================================================
% 5. HARDWARE LIMITS
% ================================================================

U1_MAX = 0.22;       % [m/s]
U2_MAX = 2.84;       % [rad/s]


%% ================================================================
% 6. SIMULATION PARAMETERS
% ================================================================

t_start = 0;
t_end   = 100;

dt = 1e-3;


%% ================================================================
% 7. INITIAL INFORMATION
% ================================================================

r20 = sum(z0.^2);

s0 = 1/(1 + c_lyap*r20);


fprintf('\n========================================\n');
fprintf(' PARAMETERS\n');
fprintf('========================================\n');

fprintf('\nGeometry:\n');
fprintf(' d0 = %.4f m\n',d0);
fprintf(' d1 = %.4f m\n',d1);
fprintf(' l1 = %.4f m\n',l1);
fprintf(' l2 = %.4f m\n',l2);

fprintf('\nNilpotent coefficients:\n');
fprintf(' alpha = %.8f\n',alpha);
fprintf(' kappa = %.8f\n',kappa);
fprintf(' theta = %.8f\n',theta);
fprintf(' beta  = %.8f\n',beta_geom);

fprintf('\nController:\n');
fprintf(' gamma   = %.6f\n',gamma);
fprintf(' c       = %.6f\n',c_lyap);
fprintf(' epsilon = %.6f s\n',epsilon);
fprintf(' omega   = %.6f rad/s\n',omega);

fprintf('\nFrequency multipliers:\n');
fprintf(' k12   = %d\n',k12);
fprintf(' k112  = %d\n',k112);
fprintf(' k1112 = %d\n',k1112);

fprintf('\nInitial physical state x0:\n');
disp(x0);

fprintf('Initial transformed state z0:\n');
disp(z0);

fprintf('||z0||       = %.8f\n',norm(z0));
fprintf('||z0||^2     = %.8f\n',r20);
fprintf('s(z0)        = %.8f\n',s0);
fprintf('effective gain gamma*s = %.8e\n',gamma*s0);


%% ================================================================
% 8. INITIAL CONTROL-BOUND DIAGNOSTIC
% ================================================================

[a1_0,a2_0,a12_0,a112_0,a1112_0] = ...
    amps( ...
        z0, ...
        alpha,kappa,theta, ...
        gamma,c_lyap);


[B1_0,B2_0, ...
 C12_0,C112_0,C1112_0] = ...
    control_bounds_from_amps( ...
        a1_0,a2_0,a12_0,a112_0,a1112_0, ...
        omega,k12,k112,k1112);


fprintf('\nInitial amplitude diagnostics:\n');

fprintf(' a1      = %+12.6e\n',a1_0);
fprintf(' a2      = %+12.6e\n',a2_0);
fprintf(' a12     = %+12.6e\n',a12_0);
fprintf(' a112    = %+12.6e\n',a112_0);
fprintf(' a1112   = %+12.6e\n',a1112_0);

fprintf('\nOscillatory contributions:\n');

fprintf(' C12         = %.6f\n',C12_0);
fprintf(' C112        = %.6f\n',C112_0);
fprintf(' 2*C1112     = %.6f\n',2*C1112_0);

fprintf('\nInitial conservative bounds:\n');

fprintf(' B1(z0) = %.6f m/s\n',B1_0);
fprintf(' B2(z0) = %.6f rad/s\n',B2_0);

fprintf('\nHardware limits:\n');

fprintf(' U1_MAX = %.6f m/s\n',U1_MAX);
fprintf(' U2_MAX = %.6f rad/s\n',U2_MAX);


%% ================================================================
% 9. RUN SIMULATION
% ================================================================

[t_sim,z_sim,u1_sim,u2_sim, ...
 B1_sim,B2_sim, ...
 t_sample,z_sample] = ...
    simulate_case( ...
        z0, ...
        alpha,kappa,theta, ...
        gamma,c_lyap, ...
        epsilon,omega, ...
        k12,k112,k1112, ...
        t_start,t_end,dt);


%% ================================================================
% 10. POST-PROCESSING
% ================================================================

norm2_sim = sum(z_sim.^2,2);

norm_sim = sqrt(norm2_sim);


% ---------------------------------------------------------------
% Actual Lyapunov function
%
% Q_c(z) = gamma/(2*c) log(1+c||z||^2)
% ---------------------------------------------------------------

Q_sim = ...
    gamma/(2*c_lyap) * ...
    log(1 + c_lyap*norm2_sim);


% Sampling-point norm
sample_norm = sqrt(sum(z_sample.^2,2));


%% ================================================================
% 11. RESULTS
% ================================================================

u1_peak = max(abs(u1_sim));

u2_peak = max(abs(u2_sim));


fprintf('\n========================================\n');
fprintf(' RESULTS\n');
fprintf('========================================\n');

fprintf('\nState norm:\n');

fprintf(' initial ||z|| = %.8f\n',norm(z0));

fprintf(' final   ||z|| = %.8f\n',norm(z_sim(end,:)));


fprintf('\nControls:\n');

fprintf(' max |u1| = %.8f m/s\n',u1_peak);

fprintf(' max |u2| = %.8f rad/s\n',u2_peak);


fprintf('\nHardware ratios:\n');

fprintf(' u1 ratio = %.4f\n',u1_peak/U1_MAX);

fprintf(' u2 ratio = %.4f\n',u2_peak/U2_MAX);


fprintf('\nConservative bounds:\n');

fprintf(' max B1 = %.8f m/s\n',max(B1_sim));

fprintf(' max B2 = %.8f rad/s\n',max(B2_sim));


if u1_peak <= U1_MAX && ...
   u2_peak <= U2_MAX

    fprintf('\n-> actual controls fit hardware limits.\n');

else

    fprintf('\n-> actual controls EXCEED hardware limits.\n');

end


%% ================================================================
% 12. PLOT STYLE
% ================================================================

col_main  = [0.00 0.20 0.70];

col_alt   = [0.85 0.33 0.10];

col_limit = [0.45 0.45 0.45];


state_labels = { ...
    'z_1', ...
    'z_2', ...
    'z_3', ...
    'z_4', ...
    'z_5' };


%% ================================================================
% FIGURE 1 -- STATE TRAJECTORIES
% ================================================================

figure( ...
    'Name','State trajectories', ...
    'NumberTitle','off');


for k = 1:5

    subplot(5,1,k);


    plot( ...
        t_sim, ...
        z_sim(:,k), ...
        'Color',col_main, ...
        'LineWidth',1.2);


    grid on;


    ylabel(state_labels{k});


    if k == 1

        title('State trajectories');

    end


    if k == 5

        xlabel('Time [s]');

    end

end


%% ================================================================
% FIGURE 2 -- CONTROL INPUTS
% ================================================================

figure( ...
    'Name','Control inputs', ...
    'NumberTitle','off');


%% ---------------------------------------------------------------
% u1
% ---------------------------------------------------------------

subplot(2,1,1);


plot( ...
    t_sim, ...
    u1_sim, ...
    'Color',col_main, ...
    'LineWidth',1.0);


hold on;


yline( ...
    U1_MAX, ...
    '--', ...
    'Color',col_limit, ...
    'LineWidth',1.2);


yline( ...
    -U1_MAX, ...
    '--', ...
    'Color',col_limit, ...
    'LineWidth',1.2);


grid on;


ylabel('u_1 [m/s]');


title('Control inputs and hardware limits');


%% ---------------------------------------------------------------
% u2
% ---------------------------------------------------------------

subplot(2,1,2);


plot( ...
    t_sim, ...
    u2_sim, ...
    'Color',col_alt, ...
    'LineWidth',1.0);


hold on;


yline( ...
    U2_MAX, ...
    '--', ...
    'Color',col_limit, ...
    'LineWidth',1.2);


yline( ...
    -U2_MAX, ...
    '--', ...
    'Color',col_limit, ...
    'LineWidth',1.2);


grid on;


ylabel('u_2 [rad/s]');

xlabel('Time [s]');


%% ================================================================
% FIGURE 3 -- CONTINUOUS AND SAMPLED STATE NORM
% ================================================================

figure( ...
    'Name','State norm', ...
    'NumberTitle','off');


plot( ...
    t_sim, ...
    norm_sim, ...
    'Color',col_main, ...
    'LineWidth',1.0);


hold on;


plot( ...
    t_sample, ...
    sample_norm, ...
    'o-', ...
    'Color',col_alt, ...
    'LineWidth',1.5, ...
    'MarkerSize',5);


grid on;


xlabel('Time [s]');

ylabel('||z(t)||');


legend( ...
    'continuous trajectory', ...
    'sampling instants', ...
    'Location','best');


title('Convergence of the state norm');


%% ================================================================
% FIGURE 4 -- SAMPLED NORM ONLY
% ================================================================

figure( ...
    'Name','Sampled state norm', ...
    'NumberTitle','off');


plot( ...
    t_sample, ...
    sample_norm, ...
    'o-', ...
    'Color',col_main, ...
    'LineWidth',1.5, ...
    'MarkerSize',6);


grid on;


xlabel('Time [s]');

ylabel('||z(j\epsilon)||');


title('State norm at sampling instants');


%% ================================================================
% FIGURE 5 -- LYAPUNOV FUNCTION
% ================================================================

figure( ...
    'Name','Lyapunov function', ...
    'NumberTitle','off');


plot( ...
    t_sim, ...
    Q_sim, ...
    'Color',col_main, ...
    'LineWidth',1.5);


grid on;


xlabel('Time [s]');

ylabel('Q_c(z)');


title(sprintf( ...
    'Q_c(z) = \\gamma/(2c) ln(1+c||z||^2),  c = %.1f', ...
    c_lyap));


%% ================================================================
% FIGURE 6 -- CONTROL ENVELOPES
% ================================================================

figure( ...
    'Name','Control envelopes', ...
    'NumberTitle','off');


subplot(2,1,1);


plot( ...
    t_sim, ...
    B1_sim, ...
    'Color',col_main, ...
    'LineWidth',1.3);


hold on;


yline( ...
    U1_MAX, ...
    '--', ...
    'Color',col_limit, ...
    'LineWidth',1.2);


grid on;


ylabel('B_1 [m/s]');

title('Conservative control envelopes');


subplot(2,1,2);


plot( ...
    t_sim, ...
    B2_sim, ...
    'Color',col_alt, ...
    'LineWidth',1.3);


hold on;


yline( ...
    U2_MAX, ...
    '--', ...
    'Color',col_limit, ...
    'LineWidth',1.2);


grid on;


ylabel('B_2 [rad/s]');

xlabel('Time [s]');


%% ================================================================
% LOCAL FUNCTIONS
% ================================================================


function g = g1_fun(z,kappa)

    g = [ ...
         1;
         0;
        -z(2);
        -z(3);
         kappa*z(2)^3-z(4) ];

end


function g = g2_fun(z,alpha,theta)

    g = [ ...
         0;
         1;
         0;
        -alpha*z(2)^2;
        -theta*z(2)*z(3) ];

end


%% ================================================================
% AMPLITUDES
%
% Q_c(z) =
%
%   gamma/(2*c) log(1+c||z||^2)
%
% grad Q =
%
%   gamma*z/(1+c||z||^2)
%
% Therefore:
%
%   gain(z) = gamma/(1+c||z||^2)
%
% and the amplitudes retain the same structure as the paper.
% ================================================================

function [a1,a2,a12,a112,a1112] = ...
    amps( ...
        z, ...
        alpha,kappa,theta, ...
        gamma,c)


    r2 = sum(z.^2);


    gain = ...
        gamma/(1 + c*r2);


    a1 = ...
        -gain*z(1);


    a2 = ...
        -gain*z(2);


    a12 = ...
        -gain * ...
        ( ...
          z(1)*z(2) ...
        + z(3) ...
        );


    a112 = ...
        -gain * ...
        ( ...
          z(1)*z(3) ...
        + z(4) ...
        + alpha*z(2)^3 ...
        );


    a1112 = ...
        -gain * ...
        ( ...
          (alpha + 2*kappa - theta) ...
             * z(1)*z(2)^3 ...
        + (alpha + 3*kappa) ...
             * z(2)^2*z(3) ...
        + z(1)*z(4) ...
        + z(5) ...
        );

end


%% ================================================================
% CONTROL u1
% ================================================================

function u1 = u1_ctrl( ...
    t,zf, ...
    omega,k12,k112,k1112, ...
    alpha,kappa,theta, ...
    gamma,c)


    [a1,~,a12,a112,a1112] = ...
        amps( ...
            zf, ...
            alpha,kappa,theta, ...
            gamma,c);


    u1 = ...
          a1 ...
        + sqrt( ...
            2*k12*omega*abs(a12)) ...
          * cos(k12*omega*t) ...
        + (2*k112*omega)^(2/3) ...
          * abs(a112)^(1/3) ...
          * cos(k112*omega*t) ...
        + 2*(k1112*omega)^(3/4) ...
          * abs(a1112)^(1/4) ...
          * cos(k1112*omega*t);

end


%% ================================================================
% CONTROL u2
% ================================================================

function u2 = u2_ctrl( ...
    t,zf, ...
    omega,k12,k112,k1112, ...
    alpha,kappa,theta, ...
    gamma,c)


    [~,a2,a12,a112,a1112] = ...
        amps( ...
            zf, ...
            alpha,kappa,theta, ...
            gamma,c);


    u2 = ...
          a2 ...
        + sqrt( ...
            2*k12*omega*abs(a12)) ...
          * sign(a12) ...
          * sin(k12*omega*t) ...
        - 2*(2*k112*omega)^(2/3) ...
          * abs(a112)^(1/3) ...
          * sign(a112) ...
          * cos(2*k112*omega*t) ...
        - 6*(k1112*omega)^(3/4) ...
          * abs(a1112)^(1/4) ...
          * sign(a1112) ...
          * sin(3*k1112*omega*t);

end


%% ================================================================
% CONTROL BOUNDS
% ================================================================

function [B1,B2,C12,C112,C1112] = ...
    control_bounds_from_amps( ...
        a1,a2,a12,a112,a1112, ...
        omega,k12,k112,k1112)


    C12 = ...
        sqrt( ...
            2*k12*omega*abs(a12));


    C112 = ...
        (2*k112*omega)^(2/3) ...
        * abs(a112)^(1/3);


    C1112 = ...
        (k1112*omega)^(3/4) ...
        * abs(a1112)^(1/4);


    B1 = ...
          abs(a1) ...
        + C12 ...
        + C112 ...
        + 2*C1112;


    B2 = ...
          abs(a2) ...
        + C12 ...
        + 2*C112 ...
        + 6*C1112;

end


%% ================================================================
% SAMPLE-AND-HOLD SIMULATION
%
% The amplitudes are computed from z(t_j) and remain frozen on
%
%     [t_j , t_j + epsilon).
%
% This is the sampling structure used in the paper.
% ================================================================

function [t_out,z_out,u1_out,u2_out, ...
          B1_out,B2_out, ...
          t_sample,z_sample] = ...
    simulate_case( ...
        z0, ...
        alpha,kappa,theta, ...
        gamma,c, ...
        epsilon,omega, ...
        k12,k112,k1112, ...
        t_start,t_end,dt)


    N = ceil( ...
        (t_end-t_start)/epsilon);


    z_cur = z0(:);


    t_out = [];

    z_out = [];


    frozen = zeros(5,N);


    B1_interval = zeros(N,1);

    B2_interval = zeros(N,1);


    t_sample = zeros(N+1,1);

    z_sample = zeros(N+1,5);


    t_sample(1) = t_start;

    z_sample(1,:) = z_cur.';


    options = odeset( ...
        'RelTol',1e-7, ...
        'AbsTol',1e-9, ...
        'MaxStep',dt);


    n_intervals = 0;


    %% ------------------------------------------------------------
    % INTEGRATION
    % ------------------------------------------------------------

    for j = 1:N


        tj = ...
            t_start + ...
            (j-1)*epsilon;


        if tj >= t_end
            break;
        end


        tj1 = ...
            min( ...
                tj+epsilon, ...
                t_end);


        % --------------------------------------------------------
        % Frozen controller state
        % --------------------------------------------------------

        zf = z_cur;


        frozen(:,j) = zf;


        n_intervals = j;


        % --------------------------------------------------------
        % Bounds for this frozen state
        % --------------------------------------------------------

        [a1,a2,a12,a112,a1112] = ...
            amps( ...
                zf, ...
                alpha,kappa,theta, ...
                gamma,c);


        [B1_interval(j), ...
         B2_interval(j)] = ...
            control_bounds_from_amps( ...
                a1,a2,a12,a112,a1112, ...
                omega,k12,k112,k1112);


        % --------------------------------------------------------
        % Closed-loop system
        % --------------------------------------------------------

        ode = @(t,z) ...
              g1_fun(z,kappa) ...
              * u1_ctrl( ...
                    t,zf, ...
                    omega,k12,k112,k1112, ...
                    alpha,kappa,theta, ...
                    gamma,c) ...
            + g2_fun(z,alpha,theta) ...
              * u2_ctrl( ...
                    t,zf, ...
                    omega,k12,k112,k1112, ...
                    alpha,kappa,theta, ...
                    gamma,c);


        [t_seg,z_seg] = ...
            ode45( ...
                ode, ...
                [tj tj1], ...
                z_cur, ...
                options);


        % --------------------------------------------------------
        % Store trajectory
        % --------------------------------------------------------

        if isempty(t_out)


            t_out = t_seg;

            z_out = z_seg;


        else


            t_out = [ ...
                t_out;
                t_seg(2:end) ];


            z_out = [ ...
                z_out;
                z_seg(2:end,:) ];


        end


        z_cur = ...
            z_seg(end,:).';


        t_sample(j+1) = ...
            tj1;


        z_sample(j+1,:) = ...
            z_cur.';


    end


    %% ------------------------------------------------------------
    % Trim
    % ------------------------------------------------------------

    frozen = ...
        frozen(:,1:n_intervals);


    B1_interval = ...
        B1_interval(1:n_intervals);


    B2_interval = ...
        B2_interval(1:n_intervals);


    t_sample = ...
        t_sample(1:n_intervals+1);


    z_sample = ...
        z_sample(1:n_intervals+1,:);


    %% ------------------------------------------------------------
    % Reconstruct controls
    % ------------------------------------------------------------

    n = numel(t_out);


    u1_out = zeros(n,1);

    u2_out = zeros(n,1);


    B1_out = zeros(n,1);

    B2_out = zeros(n,1);


    for i = 1:n


        j = ...
            floor( ...
                (t_out(i)-t_start) ...
                / epsilon) ...
            + 1;


        j = ...
            min( ...
                max(j,1), ...
                n_intervals);


        zf = frozen(:,j);


        u1_out(i) = ...
            u1_ctrl( ...
                t_out(i),zf, ...
                omega,k12,k112,k1112, ...
                alpha,kappa,theta, ...
                gamma,c);


        u2_out(i) = ...
            u2_ctrl( ...
                t_out(i),zf, ...
                omega,k12,k112,k1112, ...
                alpha,kappa,theta, ...
                gamma,c);


        B1_out(i) = ...
            B1_interval(j);


        B2_out(i) = ...
            B2_interval(j);


    end

end