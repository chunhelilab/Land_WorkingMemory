function U_min = findUmin(pca1, pca2, x, U)
U_min = zeros(1, size(x, 1));
for k = 1 : size(x, 1)
    [~, column] = min(abs(pca1(1, :)-x(k, 1)));
    [~, row] = min(abs(pca2(:, 1)-x(k, 2)));
    U_min(k) = U(row, column);
end
end