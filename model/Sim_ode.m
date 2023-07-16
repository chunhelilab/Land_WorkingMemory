function [init_value, r_end, S_end] = Sim_ode(ii, figure_folder_name, Ttotal, dt, J_min, J_max, G, varargin)
%% parameters

% default parameters
load('../model/parameters.mat');

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


%% initialization

% Number of time points
NT = floor(Ttotal / dt);

% initialization
S_A = rand(area_num, 1);
S_B = rand(area_num, 1);
S_C = rand(area_num, 1);

init_value = [S_A; S_B; S_C];

t = (1 : NT)' .* dt;
r_A_record = zeros(NT, area_num);
r_B_record = zeros(NT, area_num);
r_C_record = zeros(NT, area_num);

S_A_record = zeros(NT, area_num);
S_B_record = zeros(NT, area_num);
S_C_record = zeros(NT, area_num);


%% simulation
for i_t = 1 : NT
    
    I_A_net = WFF * S_A;
    I_B_net = WFF * S_B;
    I_C_net = WFB * (S_A + S_B);
    
    I_A = J_s .* S_A + J_c .* S_B + J_EI .* S_C + I_0A + I_ramp;
    I_B = J_c .* S_A + J_s .* S_B + J_EI .* S_C + I_0B + I_ramp;
    I_C = J_IE .* S_A + J_IE .* S_B + J_II .* S_C + I_0C + I_ramp;
    
    I_A = I_A + I_A_net;
    I_B = I_B + I_B_net;
    I_C = I_C + I_C_net;
    
    I_A(1, 1) = I_A(1, 1) + I_stim;  % visual stimulus to the population A of V1 area
    I_B(1, 1) = I_B(1, 1) + I_dist;  % visual stimulus to the population B of V1 area
    I_C(spine_count_idx(27:30), 1) = I_C(spine_count_idx(27:30), 1) + I_inhi; % inhibitory current to the population C of F7,8B, 9/46v, 9/46d areas
    
    r_A = (a .* I_A - b) ./ (1 - exp(-d .* (a .* I_A - b)));
    r_B = (a .* I_B - b) ./ (1 - exp(-d .* (a .* I_B - b)));
    r_C = (c_1 .* I_C - c_0) ./ g_I + r_0;
    r_C(r_C < 0) = 0;
    
    S_A = S_A + (-S_A ./ tau_N + gamma_E .* (1 - S_A) .* r_A) * dt;
    S_B = S_B + (-S_B ./ tau_N + gamma_E .* (1 - S_B) .* r_B) * dt;
    S_C = S_C + (-S_C ./ tau_G + gamma_I .* r_C) .* dt;
    
    
    S_A_record(i_t, :) = S_A';
    S_B_record(i_t, :) = S_B';
    S_C_record(i_t, :) = S_C';
    
    r_A_record(i_t, :) = r_A';
    r_B_record(i_t, :) = r_B';
    r_C_record(i_t, :) = r_C';
    
end
r_end = [r_A_record(end, :), r_B_record(end, :), r_C_record(end, :)];
S_end = [S_A_record(end, :), S_B_record(end, :), S_C_record(end, :)];

%% plot
h = figure(ii);
set(h, 'position', [10 10 950 700], 'Visible', 'off');
start_time = dt / dt;
end_time = NT;
step = 100;
for k = 1 : area_num
    subplot(5, 6, k)
    plot(t(start_time:step:end_time), r_A_record(start_time:step:end_time, spine_count_idx(k)), 'Color', '#4D21BD', 'LineWidth', 1.5)
    hold on
    plot(t(start_time:step:end_time), r_B_record(start_time:step:end_time, spine_count_idx(k)), 'Color', '#59A95A', 'LineWidth', 1.5)
    %         hold on
    %         plot(t(start_time:step:end_time), r_C_record(start_time:step:end_time, spine_count_idx(k)), 'Color', '#FF903D', 'LineWidth', 1)
    title(area_name{k}, 'FontSize', 12)
    if k == 1
        legend('r_A', 'r_B', 'box', 'off')
    end
    if mod(k, 6) == 1
        ylabel('Rate (sp/s)')
    end
    if k > 24
        xlabel('Time (s)')
    end
    ylim([0, 60])
end
print(h, [figure_folder_name, '\figure', int2str(ii)], '-dpng', '-r180');
close(h)

end





