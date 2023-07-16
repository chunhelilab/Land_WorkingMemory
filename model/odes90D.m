function F = odes90D(J_min, J_max, G, syms_90D, varargin)
%%
% default parameters
load('parameters.mat');

p = inputParser;            % 函数的输入解析器
addParameter(p, 'J_c', 0.0107);  % 设置变量名和默认参数
addParameter(p, 'J_EI', -0.31);      
addParameter(p, 'J_II', -0.12);
addParameter(p, 'I_stim', 0);
addParameter(p, 'I_dist', 0);
addParameter(p, 'I_ramp', 0);
addParameter(p, 'Perturbation_area', []);
addParameter(p, 'I_inhi', 0);
addParameter(p, 'LWM', false);
parse(p, varargin{:});       % 对输入变量进行解析，如果检测到前面的变量被赋值，则更新变量取值

J_c =  p.Results.J_c;
J_EI = p.Results.J_EI;
J_II = p.Results.J_II;
I_stim = p.Results.I_stim;
I_dist = p.Results.I_dist;
I_ramp = p.Results.I_ramp;
Perturbation_area = p.Results.Perturbation_area;
I_inhi = p.Results.I_inhi;
LWM = p.Results.LWM;

% localized working memory
if LWM == true
    load('../Processed_Data/FLN_matrix.mat', 'FLN');  % (to, from)
    FLN = tril(FLN); % for localized working memory
    FLN = k_1 .* FLN .^ (k_2);
    norm = sum(FLN, 2);
    for k = 1 : area_num  % for localized working memory
        FLN(k,:) = FLN(k,:) / max(norm(k, 1), 0.1);
    end
end

for k = Perturbation_area
    FLN(:, k) = 0; % perturbation
    FLN(k, :) = 0; % perturbation
end

% other parameters
zeta = tau_G * gamma_I * c_1 / (g_I - J_II * tau_G * gamma_I * c_1);
J_0 = J_s + J_c + 2 * J_EI * J_IE * zeta;  % (nA)

J_s = J_min + (J_max - J_min) * spine_count;  % (nA)
J_IE = (J_0 - J_s - J_c) ./ (2 * J_EI * zeta);  % (nA)
Z = 2 * c_1 * tau_G * gamma_I * J_EI / (c_1 * tau_G * gamma_I * J_II - g_I);
WFF = zeros(area_num, area_num);
WFB = zeros(area_num, area_num);
for k = 1 : area_num
    WFF(k, :) = J_s(k) ./ max(J_s) .* FLN(k, :) .* G;
    WFB(k, :) = J_IE(k) ./ max(J_IE) .* FLN(k, :) .* G ./ Z;
end
WFF = WFF .* SLN;
WFB = WFB .* SLN_i;

%% variables
syms_S_A = syms_90D(1:area_num, :);
syms_S_B = syms_90D(area_num+1:2*area_num, :);
syms_S_C = syms_90D(2*area_num+1:3*area_num, :);


%% ODEs
I_A_net = WFF * syms_S_A;
I_B_net = WFF * syms_S_B;
I_C_net = WFB * (syms_S_A + syms_S_B);

I_A = J_s .* syms_S_A + J_c .* syms_S_B + J_EI .* syms_S_C + I_0A + I_A_net + I_ramp;
I_B = J_c .* syms_S_A + J_s .* syms_S_B + J_EI .* syms_S_C + I_0B + I_B_net + I_ramp;
I_C = J_IE .* syms_S_A + J_IE .* syms_S_B + J_II .* syms_S_C + I_0C + I_C_net + I_ramp;

I_A(1) = I_A(1) + I_stim; % V1, stimulus
I_B(1) = I_B(1) + I_dist; % V1, distractor
I_C(spine_count_idx(27:30), 1) = I_C(spine_count_idx(27:30), 1) + I_inhi; % F7,8B, 9/46v, 9/46d, inhibitory current

r_A = (a .* I_A - b) ./ (1 - exp(-d .* (a .* I_A - b)));
r_B = (a .* I_B - b) ./ (1 - exp(-d .* (a .* I_B - b)));
r_C = (c_1 .* I_C - c_0) ./ g_I + r_0;


%% dirving force F
F = [-syms_S_A ./ tau_N + gamma_E .* (1 - syms_S_A) .* r_A; -syms_S_B ./ tau_N + gamma_E .* (1 - syms_S_B) .* r_B; -syms_S_C ./ tau_G + gamma_I .* r_C];

end