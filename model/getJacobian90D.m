function Jacob_F = getJacobian90D(x)
%% Jacobian matrix
mf_Jacob_F_90D = matlabFunction(jacobian(F, syms_90D));
load mf_Jacob_F_90D

Jacob_F = zeros(size(Jacob_F_sys, 1), size(Jacob_F_sys, 2), size(x, 2));
for k = 1 : size(x, 2)
    Jacob_F(:, :, k) = mf_Jacob_F_90D(x(1, k), x(2, k), x(3, k), x(4, k), x(5, k), x(6, k), x(7, k), x(8, k), x(9, k), x(10, k), x(11, k), x(12, k), x(13, k), x(14, k), x(15, k), x(16, k), x(17, k), x(18, k), x(19, k),...
    x(20, k), x(21, k), x(22, k), x(23, k), x(24, k), x(25, k), x(26, k), x(27, k), x(28, k), x(29, k),x(30, k), x(31, k), x(32, k), x(33, k), x(34, k), x(35, k), x(36, k), x(37, k), x(38, k), x(39, k), ...
    x(40, k), x(41, k), x(42, k), x(43, k), x(44, k), x(45, k), x(46, k), x(47, k), x(48, k), x(49, k), x(50, k), x(51, k), x(52, k), x(53, k), x(54, k), x(55, k), x(56, k), x(57, k), x(58, k), x(59, k), ...
    x(60, k), x(61, k), x(62, k), x(63, k), x(64, k), x(65, k), x(66, k), x(67, k), x(68, k), x(69, k), x(70, k), x(71, k), x(72, k), x(73, k), x(74, k), x(75, k), x(76, k), x(77, k), x(78, k), x(79, k), ...
    x(80, k), x(81, k), x(82, k), x(83, k), x(84, k), x(85, k), x(86, k), x(87, k), x(88, k), x(89, k), x(90, k));
end

end
