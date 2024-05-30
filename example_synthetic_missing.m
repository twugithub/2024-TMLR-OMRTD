clearvars;
close all;
clc

rng('shuffle');
addpath('synthetic');
addpath('OMRTD');

n1 = 50;
num_samples = 2000;
n3 = 20;
rho = 0.1;

lambda1 = 1/sqrt(n1);
lambda2 = 1/sqrt(n1);

rank_ratio = 0.2;
d = round(n1 * rank_ratio);
[ U, V, Z, E ] = gen_tensor_data( n1, n3, num_samples, rank_ratio, rho );
Z = Z + E;

eta = 0.8;
Z_miss = NaN(size(Z));

for i = 1:size(Z,2)
    tempdata = squeeze(Z(:,i,:));
    chosen = randperm(size(tempdata,1)*size(tempdata,2), round(eta*size(tempdata,1)*size(tempdata,2)));
    tmp = NaN(size(tempdata,1),size(tempdata,2));
    tmp(chosen) = tempdata(chosen);
    Z_miss(:,i,:) = reshape(tmp, [size(tempdata,1),1,size(tempdata,2)]);
end

tic

[L_ROMRTD, R_ROMRTD, E_ROMRTD] = solve_romrtd( Z_miss, d, lambda1, lambda2 );
Ln = orth_tensor(L_ROMRTD{num_samples+1});
ev = compute_EV(Ln, U);

toc