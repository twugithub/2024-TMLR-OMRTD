function [L_est, R_est, E_est] = solve_omrtd( Z, d, lambda1, lambda2 )

%% initialization
[n1,n2,n3] = size(Z);
L_est = cell(1, n2+1);
L_est{1} = rand(n1,d,n3);
R_est = zeros(n2,d,n3);
E_est = zeros(n1,n2,n3);

A = zeros(d,d,n3);
B = zeros(n1,d,n3);

%% online optimization
for t = 1:n2

    if mod(t, 200) == 0
        fprintf('OMRTD: Access sample %d\n', t);
    end

    z = Z(:,t,:);

    [r, e] = solve_re(z, L_est{t}, lambda2);

    R_est(t,:,:) = tran(r);
    E_est(:,t,:) = e;

    A = A + tprod(r, tran(r));
    B = B + tprod(z-e, tran(r));
    L_est{t+1} = solve_L(L_est{t}, A, B, lambda1);

end

end