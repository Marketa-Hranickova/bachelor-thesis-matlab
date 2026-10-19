function ok = satisfies_c1(k12, k112, k1112)
%SATISFIES_C1  Condition C1 of Theorem 1 (Zuyev & Grushkovskaya, 2403.03341).
%
%   C1) there are no resonances of order up to 4 in each of the tuples
%       (k12, k112, k1112), (k12, 2*k112, k1112),
%       (k12, k112, 3*k1112), (k12, 2*k112, 3*k1112)
%
%   A tuple (k_1,...,k_s) is in resonance of order K if there exist SETWISE
%   COPRIME integers (c_1,...,c_s) with sum|c_i| = K and sum c_i*k_i = 0.
%
%   Smaller multipliers mean smaller commanded velocities, so it is tempting
%   to shrink them -- this is the condition that says how far you may go.

    tuples = [k12, k112,   k1112;
              k12, 2*k112, k1112;
              k12, k112,   3*k1112;
              k12, 2*k112, 3*k1112];

    ok = true;
    for r = 1:size(tuples, 1)
        if has_resonance(tuples(r, :), 4)
            ok = false;
            return
        end
    end
end


function res = has_resonance(kt, max_order)
%HAS_RESONANCE  True if the tuple has a resonance of order <= max_order.
    res = false;
    for c1 = -max_order:max_order
        for c2 = -max_order:max_order
            for c3 = -max_order:max_order
                K = abs(c1) + abs(c2) + abs(c3);
                if K == 0 || K > max_order
                    continue
                end
                if gcd(gcd(abs(c1), abs(c2)), abs(c3)) ~= 1
                    continue        % not setwise coprime
                end
                if c1*kt(1) + c2*kt(2) + c3*kt(3) == 0
                    res = true;
                    return
                end
            end
        end
    end
end
