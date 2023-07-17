clear; clc;
addpath(genpath('./'))
tic;
%% Initialize Matlab Parallel Computing Enviornment
init_parallel(6)

%% calculate the mean and variance of potential landscape

load('./model/parameters.mat')
file_name = 'data_G0.48_Js0.21-0.3_I0.00_Stim.mat';
load(['./simulation/data/', file_name])

area_num = 30;
J_min = str2num(file_name(strfind(file_name, 's')+1 : (strfind(file_name, '-')-1)));  % (nA) the minimum recurrent connectivity strength
J_max = str2num(file_name(strfind(file_name, '-')+1 : (strfind(file_name, 'I')-2)));  % (nA) the maximum recurrent connectivity strength
G = str2num(file_name(strfind(file_name, 'G')+1 : (strfind(file_name, 'J')-2)));  % global coupling factor
J_II=-0.12; I_dist = 0; I_ramp = 0; Perturbation_area = []; I_inhi = 0; LWM = false;
I_stim = str2num(file_name(strfind(file_name, 'I')+1 : strfind(file_name, '_Stim')-1));


%%%%% determine the unique attractors
var = S_end_all;
[uni_var, ~, IC] = uniquetol(var, 0.1, 'ByRows', true);
weight = hist(IC, unique(IC)) /  size(var, 1);
[weight, I] = sort(weight, 'descend');
uni_var = uni_var(I, :);

mu = uni_var;
nAttractor = length(weight);



%%%%% Jacobian matrix of ODEs
sig = cell(1, nAttractor);
for k = 1 : area_num
    syms_90D(k, 1) = sym (['S_A_', num2str(k)]);
    syms_90D(area_num+k, 1) = sym (['S_B_', num2str(k)]);
    syms_90D(2*area_num+k, 1) = sym (['S_C_', num2str(k)]);
end
F = odes90D(J_min, J_max, G, syms_90D, 'J_II', J_II, 'I_stim', I_stim, 'I_dist', I_dist, 'I_ramp', I_ramp, 'Perturbation_area', Perturbation_area, 'I_inhi', I_inhi, 'LWM', LWM);
Jacob_F_sys = jacobian(F, syms_90D);


%%%%% variance matrix at steady state
B = -2 * diag([repmat(sigma_A/tau_noise, 1, area_num), repmat(sigma_B/tau_noise, 1, area_num), repmat(sigma_C/tau_noise, 1, area_num)]);

for k = 1 : nAttractor
    Jacob_F = subs(Jacob_F_sys, syms_90D, mu(k, :)');
    Jacob_F = double(Jacob_F);
    sig{1, k} = calculateVariance(Jacob_F, B);
    if find(diag(sig{k}) < 0)
        error('The variance is negative!')
    end
end

%% dimension reduction
[eigVector, eigValue] = dimensionReduction(mu, weight, sig);
disp(['The first two componant explains ', num2str(sum(eigValue(1:2, 1))/sum(eigValue)*100), '% variance'])

mu_pca = zeros(nAttractor, 2);
sig_pca = cell(nAttractor, 1);
for k = 1 : nAttractor
    mu_pca(k, :) = eigVector(:, 1:2)' * mu(k, :)';
    sig_pca{k} = eigVector(:, 1:2)' * sig{k} * eigVector(:, 1:2);
end

% plot eigenvalue
figure(1)
semilogx(eigValue, '-*', 'LineWidth', 1.5)
ylabel('Eigenvalue')

% plot eigenvector
figure(2);
set(gcf,'outerposition', [100 100 1000 650]);
for k = 1 : 2
    subplot(2, 1, k)
    b = bar(reshape(eigVector(:, k), [30, 3]), 'EdgeColor', [1, 1, 1]);
    b(1).FaceColor = '#4D21BD';
    b(2).FaceColor = '#59A95A';
    b(3).FaceColor = '#FF903D';
    set(gca,'FontSize',8)
    box off
    xlabel('Variables','FontSize', 10)
    xticks([1:1:30])
    xticklabels(area_name)
    xtickangle(60)
    ylabel('A.U.','FontSize', 10)
    legend('S_A', 'S_B', 'S_C', 'box', 'off', 'Location', 'northwest')
    title(['Eigenvector ', num2str(k)],'FontSize', 12)
end

% Rerank attractors according to ascending PC1
[mu_pca, C] = sortrows(mu_pca);
sig_pca = sig_pca(C);
mu = mu(C, :);
sig = sig(C);
weight = weight(C);
disp(mu_pca)

% calculate the landscape
y_max = [4, 4]; %% Range of the landscape
y_min = [-4, -1];
g = 501;
step = (y_max - y_min) / (g - 1); %% Length of the step
[pca1, pca2] = meshgrid(y_min(1) : step(1) : y_max(1), y_min(2) : step(2) : y_max(2)); %% Grid
P_ss = zeros(g);
for k = 1 : size(mu_pca, 1)
    for m = 1 : g
        for n = 1 : g
            P_ss(m, n) =  P_ss(m, n) + weight(k) * multivariate_normal_distribution([pca1(m, n); pca2(m, n)], mu_pca(k, :)', sig_pca{k}, 2);
        end
    end
end
P_ss = P_ss / sum(sum(P_ss)); % normalization
U = -log(P_ss);
U(U > 100)  = 100;

%% Minimum Action Path
Func = @(syms_90D)odes90D(J_min, J_max, G, syms_90D, 'J_II', J_II, 'I_stim', I_stim, 'I_dist', I_dist, 'I_ramp', I_ramp, 'Perturbation_area', Perturbation_area, 'I_inhi', I_inhi, 'LWM', LWM);  % Force
params.Dimension = size(mu, 2); % The dimension of the system.
dFunc = []; % Jaconbian matrix

num_path = nAttractor*(nAttractor-1);
action= zeros(num_path, 1);
action_record = cell(num_path, 1);
ycell = cell(num_path, 1);
start_attractor = repelem([1:1:nAttractor], nAttractor-1);
end_attractor = [];
for i = 1 : 1 : nAttractor
    temp = 1 : 1 : nAttractor;
    temp(temp == i) = [];
    end_attractor = [end_attractor, temp];
end

parfor i = 1 : num_path  % this parfor loop requires about 2.5h to finish when using 6 cores
    Spot1 = start_attractor(i); % Initial point
    Spot2 = end_attractor(i); % End point
    init = mu([Spot1, Spot2], :)';
    [ActionVal, Path] = minActionPath(params, init, Func, dFunc);
    action_record{i} = ActionVal;
    action(i) = ActionVal(end, 1);
    ycell(i) = {[mu(Spot1, :)' Path mu(Spot2, :)']};
end

% Convergence of action
figure(3)
for i = 1 : nAttractor
    for j = 1 : nAttractor
        if i ~= j
            subplot(nAttractor, nAttractor, (i-1)*nAttractor+j)
            temp = find((start_attractor == i) & (end_attractor == j));
            plot(1:size(action_record{temp}, 1), action_record{temp}(1:size(action_record{temp},1), :), 'LineWidth', 2);
            xlabel('Iteration Number')
            ylabel('Action')
        end
    end
end

% plot discrete path
figure(4);
population_idx = repelem({' A', ' B', ' C'}, area_num);
path_temp = ycell{2}; % MA to MB
image(path_temp', 'CDataMapping', 'scaled')
box off
load('./utils/coolwarm_cmap.mat')
colormap(coolwarm_cmap)
x_tick_label = [];
for k = 1 : 90
    x_tick_label = [x_tick_label, strcat(area_name{k-floor((k-1)/area_num)*area_num}, population_idx(k))];
end
set(gca, 'XTick', [1 : 90], 'XTickLabel',  x_tick_label, 'YTick', [1 : 6 : params.N+1], 'YTickLabel', [0 : params.TMax/params.N*6 : params.TMax], 'fontsize', 7);
xtickangle(90)
xlabel('Brain areas', 'fontsize', 14);
ylabel('Time', 'fontsize', 14);
colorbar();


% Calculate the paths after dimension reduction
path_pca = cell(num_path, 1);
for i = 1 : num_path
    path_pca{i} =  eigVector(:, 1:2)'*ycell{i};
end

%% Plotting

% Plot the landscape
figure(5)
clf;
surf(pca1, pca2, U);
shading interp
xlabel('PC1', 'FontSize', 14)
ylabel('PC2', 'FontSize', 14)
zlabel('U', 'FontSize', 14)
axis([y_min(1) y_max(1) y_min(2) y_max(2) 0 Inf])
savetitle = ['G=', num2str(G), ' J_{min}=', num2str(J_min), ' J_{max}=', num2str(J_max), ' I_{Stim}=', num2str(I_stim)];
title(savetitle, 'FontSize', 16)
view(40, 80)
hold on
caxis([0 100]);

% Plot the minima of the potential 
U_min = findUmin(pca1, pca2, mu_pca, U);
for k = 1 : size(mu_pca, 1)
    hold on
    scatter3(mu_pca(k, 1), mu_pca(k, 2), U_min(k)+2, 40, 's', 'filled', 'MarkerEdgeColor', '#FB9A99', 'MarkerFaceColor', '#FB9A99')
end

% Plot the grid
for i = 1 : floor(size(pca1,1)/20)
    plot3(pca1(20*i-1,:), pca2(20*i-1,:), U(20*i-1,:), 'Color', [0.4 0.4 0.4], 'LineWidth', 0.2);
end
for i = 1 : floor(size(pca1,2)/20)
    plot3(pca1(:,20*i-1), pca2(:,20*i-1), U(:,20*i-1), 'Color', [0.4 0.4 0.4], 'LineWidth',0.2);
end

% Plot the paths
for i = 2
    for j = 1
        temp = find((start_attractor == i) & (end_attractor == j));
        z3path = griddata(pca1, pca2, U, path_pca{temp}(1,:), path_pca{temp}(2,:));
        plot3(path_pca{temp}(1,:), path_pca{temp}(2,:), z3path+7, 'w', 'LineWidth', 2);
        
        temp = find((start_attractor == j) & (end_attractor == i));
        z3path = griddata(pca1, pca2, U, path_pca{temp}(1,:), path_pca{temp}(2,:));
        plot3(path_pca{temp}(1,:), path_pca{temp}(2,:), z3path+7, 'Color', [0.85,0.43,0.83], 'LineWidth',2);
    end
end
view([-25 75])
set(gcf,'outerposition', [100 100 800 650]);

% find minimum and saddle point of U
if size(mu_pca, 1) == 2
    [U_saddle, saddle_point] = findUsaddle(pca1, pca2, U, U_min);
elseif size(mu_pca, 1) == 3
    U1 = U;
    U1(:, (g+1)/2:g) = 100;
    U2 = U;
    U2(:, 1:(g-1)/2) = 100;
    [U_saddle1, saddle_point1] = findUsaddle(pca1, pca2, U1, U_min(1:2));
    [U_saddle2, saddle_point2] = findUsaddle(pca1, pca2, U2, U_min(2:3));
    U_saddle = [U_saddle1, U_saddle2];
    saddle_point = [saddle_point1; saddle_point2];
end

% Plot the saddle point
for k = 1 : size(saddle_point, 1)
    hold on
    scatter3(saddle_point(k, 1), saddle_point(k, 2), U_saddle(k)+1, 40, 'o', 'filled', 'MarkerEdgeColor', '#00FFFF', 'MarkerFaceColor', '#00FFFF')
end

%%
toc;