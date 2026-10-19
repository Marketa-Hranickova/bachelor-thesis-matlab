function [h, name] = lyapunov_h(kind, par)
%LYAPUNOV_H  State-dependent gain factor h(z) of a radial Lyapunov function.
%
%   [h, name] = lyapunov_h('quadratic')
%   [h, name] = lyapunov_h('log',  c)
%   [h, name] = lyapunov_h('sqrt', delta)
%
%   For V(z) = phi(r), r = ||z||, the controller amplitudes are
%       a(z) = -gamma * h(z) * G^-1(z) * z ,      h(z) = phi'(r)/(gamma*r).
%
%     kind         V(z)                              h(z)
%     'quadratic'  gamma*r^2/2                       1
%     'log'        gamma/(2c) * log(1 + c*r^2)       1/(1 + c*r^2)
%     'sqrt'       gamma*(sqrt(r^2 + d^2) - d)       1/sqrt(r^2 + d^2)
%
%   h is returned as a function handle  h(z)  taking the 5x1 state z.
%   'name' is a short label for printouts and plot legends.
%
%   Intuition: h(z) is how hard the controller pushes at distance r.
%   'log' pushes weakly far from the goal (small controls) and fully
%   (h -> 1) near it; 'sqrt' pushes MORE and more as r -> 0 (up to 1/delta).

    if nargin < 2, par = []; end

    switch kind
        case 'quadratic'
            h    = @(z) 1;
            name = 'quadratic';

        case 'log'
            c    = par;
            h    = @(z) 1/(1 + c*norm(z)^2);
            name = sprintf('log, c = %g', c);

        case 'sqrt'
            d    = par;
            h    = @(z) 1/sqrt(norm(z)^2 + d^2);
            name = sprintf('sqrt, delta = %g', d);

        otherwise
            error('lyapunov_h:unknownKind', ...
                  'Unknown kind ''%s'' (use quadratic | log | sqrt).', kind);
    end
end
