function [m, r, e] = solve_missing_re(z, L, lambda2)

[n1,d,n3] = size(L);

r = zeros(d,1,n3);
e = zeros(n1,1,n3);
y = zeros(n1,1,n3);
oldr = r;olde = e;

converged = false;
mu = 0.1;
rho = 1.9;
mu_max = 1e10;
maxIter = 100;
iter = 0;

I = eye(d,d);
Lfft = fft(L,[],3);
LLtinv = zeros(d,n1,n3);
for i = 1:n3
    LLtinv(:,:,i) = (Lfft(:,:,i)' * Lfft(:,:,i) + 0.01 * I) \ Lfft(:,:,i)';
end

zsq = squeeze(z);
Indicator = ~isnan(zsq(:));
missloc = isnan(zsq(:));

while ~converged
    iter = iter + 1;

    zysq = squeeze(z - y/mu);
    msq = zeros(size(zsq));
    tmp = squeeze(tprod(L, r) + e);
    msq(Indicator) = (tmp(Indicator) + mu*zysq(Indicator)) / (mu + 1);
    msq(missloc) = tmp(missloc);

    m = reshape(msq, [size(z,1), 1, size(z,3)]);

    diff_me = fft(m - e,[],3);
    r0fft = zeros(d,1,n3);
    for i = 1:n3
        r0fft(:,:,i) = LLtinv(:,:,i) * diff_me(:,:,i);
    end
    norm_r0 = sqrt(squeeze(sum(conj(r0fft).*r0fft,1)));

    for i = 1:n3
        if norm_r0(i) <= 1
            r(:,:,i) = r0fft(:,:,i);
        else
            [r(:,:,i), ~] = solve_r(squeeze(Lfft(:,:,i)), squeeze(diff_me(:,:,i)));
        end
    end
    r = ifft(r,[],3);
    e = sign(m - tprod(L, r)) .* max(abs(m - tprod(L, r)) - lambda2, 0);
    dd = squeeze(m - z);
    dd(missloc) = 0;

    stopc = max([sqrt(sum(sum(dd.^2))), sqrt(sum(sum((r - oldr).^2))), sqrt(sum(sum((e - olde).^2)))]) / (n1*n3);

    if stopc < 1e-6 || iter > maxIter
        converged = true;
    else
        ysq = squeeze(y);
        msq = squeeze(m);
        ysq(Indicator) = ysq(Indicator) + mu*(msq(Indicator) - zsq(Indicator));
        y = reshape(ysq, [size(z,1), 1, size(z,3)]);
        mu = min(rho*mu, mu_max);
        oldr = r;olde = e;
    end

end

end