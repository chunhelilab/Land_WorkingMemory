function sig = calculateVariance(Jacob_F, B)
%% The variance of Gaussian distribution, A*sigma + sigma*A' + 2*d*I
N = length(Jacob_F); % the dimension of the system
B = reshape(B, N^2, 1);
L = zeros(N^2, N^2);
for k = 1 : N
    L(((k-1)*N+1): k*N, ((k-1)*N+1): k*N) = Jacob_F;
end
for m = 1 : N
    for n = 1 : N
        L(((m-1)*N+1) : m*N, ((n-1)*N+1) : n*N) = L(((m-1)*N+1) : m*N, ((n-1)*N+1) : n*N) + Jacob_F(m, n) * eye(N);
    end
end
sig = L \ B;
sig = reshape(sig, N, N)';

end
