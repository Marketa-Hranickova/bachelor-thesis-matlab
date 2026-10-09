%UNPACK_PARAMS  Convenience for the COMMAND WINDOW only.
%
%   >> P = params();  unpack_params
%
%   Do NOT call this from inside a function. A script called from a function
%   body assigns its variables at run time, which the parser never sees, so
%   a later bare `alpha` is still compiled as the built-in. Inside a function,
%   write the assignments out in the function itself.

d0 = P.d0;  d1 = P.d1;  l1 = P.l1;  l2 = P.l2;
alpha = P.alpha;  kappa = P.kappa;  theta = P.theta;  beta = P.beta;
gamma = P.gamma;  epsilon = P.epsilon;  omega = P.omega;
k12 = P.k12;  k112 = P.k112;  k1112 = P.k1112;
z0 = P.z0;  t_start = P.t_start;  t_end = P.t_end;  dt = P.dt;
U1_MAX = P.u1_max;  U2_MAX = P.u2_max;
