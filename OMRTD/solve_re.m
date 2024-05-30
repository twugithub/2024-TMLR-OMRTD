function [r, e] = solve_re(z, L, lambda2)

[n1,d,n3] = size(L);

r = zeros(d,1,n3);
e = zeros(n1,1,n3);

converged = false;
maxIter = 100;
iter = 0;

I = eye(d,d);
Lfft = fft(L,[],3);
LLtinv = zeros(d,n1,n3);
for i = 1:n3
    LLtinv(:,:,i) = (Lfft(:,:,i)' * Lfft(:,:,i) + 0.01 * I) \ Lfft(:,:,i)';
end

while ~converged
    iter = iter + 1;

    orgr = r;
    orge = e;

    diff_ze = fft(z - e,[],3);
    r0fft = zeros(d,1,n3);
    for i = 1:n3
        r0fft(:,:,i) = LLtinv(:,:,i) * diff_ze(:,:,i);
    end
    norm_r0 = sqrt(squeeze(sum(conj(r0fft).*r0fft,1)));

    for i = 1:n3
        if norm_r0(i) <= 1
            r(:,:,i) = r0fft(:,:,i);
        else
            [r(:,:,i), eta] = solve_r(squeeze(Lfft(:,:,i)), squeeze(diff_ze(:,:,i)));
        end
    end
    r = ifft(r,[],3);
    e = sign(z - tprod(L, r)) .* max(abs(z - tprod(L, r)) - lambda2, 0);

    stopc = max(sqrt(sum(sum((r - orgr).^2))), sqrt(sum(sum((e - orge).^2)))) / (n1*n3);

    if stopc < 1e-6 || iter > maxIter
        converged = true;
    end
end

end