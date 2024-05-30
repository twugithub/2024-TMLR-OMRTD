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

tic

[L_OMRTD, R_OMRTD, E_OMRTD] = solve_omrtd( Z, d, lambda1, lambda2 );
Ln = orth_tensor(L_OMRTD{num_samples+1});
ev = compute_EV(Ln, U);

toc