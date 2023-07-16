function [V, D] = dimensionReduction(mu, weight, sig)
attractor_num = length(weight);

%% mean for multistability
Mu = weight * mu;  

%% variance for multistability
Sigma = -Mu' * Mu;
for k = 1 : attractor_num
    Sigma = Sigma + weight(k) * (sig{k} + mu(k, :)' * mu(k, :));
end

%% the first two componant
[V, D] = eig(Sigma);
[D, ind] = sort(diag(D), 'descend');
V = V(:, ind);

% make sure most of the elements of the eigenvector are positive.
for k = 1 : 2
    if sum(V(:, k)) < 0
        V(:, k) = -V(:, k);
    end
end

end
