function [r, eta] = solve_r(L, diff_ze)

aux_data = L' * diff_ze;
LtL = L' * L;

d = size(L,2);
I = eye(d,d);

start_pt = 0;
end_pt = 10;
eta = end_pt;
converged = false;

% find a correct region
while ~converged
    r = (LtL + eta * I) \ aux_data;
    nrm = norm(r);

    if nrm < 1
        converged = true;
    else
        eta = eta * 2;
    end
end

end_pt = eta;
converged = false;
eta = (start_pt + end_pt) /2;

while ~converged
    r = (LtL + eta * I) \ aux_data;
    nrm = norm(r);

    mid_pt = (start_pt + end_pt) /2;

    if abs(1 - nrm) < 1e-3
        converged = true;
    else
        if nrm > 1
            start_pt = mid_pt;
        else
            end_pt = mid_pt;
        end
        eta = (start_pt + end_pt) / 2;
    end
end

end