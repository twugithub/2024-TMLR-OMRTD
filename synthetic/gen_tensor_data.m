function [ U, V, Z, E ] = gen_tensor_data( n1, n3, num_samples, rank_ratio, rho )

d = round(n1 * rank_ratio);

U = randn(n1,d,n3);
U = orth_tensor(U);
V = randn(num_samples,d,n3);

Z = tprod(U, tran(V));

num_elements = n1 * num_samples * n3;
temp = randperm(num_elements);
numCorruptedEntries = round(rho * num_elements);
corruptedPositions = temp(1:numCorruptedEntries);
E = zeros(n1, num_samples, n3);
E(corruptedPositions) = 20 *(rand(numCorruptedEntries, 1) - 0.5);

end