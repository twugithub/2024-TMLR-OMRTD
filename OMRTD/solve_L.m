function [L_est] = solve_L(L, A, B, lambda1)

[~,d,n3] = size(L);

Afft = fft(A,[],3);
Bfft = fft(B,[],3);
Lfft = fft(L,[],3);

U = subgrad_2_inf(Lfft);

for k = 1:n3
    for j = 1:d
        bj = Bfft(:,j,k);
        lj = Lfft(:,j,k);
        aj = Afft(:,j,k);
        uj = U(:,j,k);
        temp = lj - ((Lfft(:,:,k)*aj - bj) / n3 + lambda1 * uj) / (Afft(j,j,k) + 0.001);
        Lfft(:,j,k) = temp/max(norm(temp),1);
    end
end

L_est = ifft(Lfft,[],3);

end



function [U] = subgrad_2_inf(L)

[n1,d,n3] = size(L);
U = zeros(n1,d,n3);
norm2 = zeros(n1,n3);

for i = 1:n1
    for k = 1:n3
        norm2(i,k) = norm(L(i,:,k));
    end
end

[Q1, Q2] = find(norm2 == max(norm2(:)));

for k = 1:n3
    mu = zeros(n1,1);
    idx = find(Q2 == k);
    Q = Q1(idx);
    mu(Q) = 1 / length(Q1);
    mu = repmat(mu, 1, d);
    U(:,:,k) = mu .* L(:,:,k);
end

end