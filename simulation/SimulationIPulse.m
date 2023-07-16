clear;
clc;
tic;
addpath(genpath('../'))
%% Initialize Matlab Parallel Computing Enviornment
init_parallel(6)

%% area name and number
load('../model/parameters.mat', 'area_name', 'area_num');

%%
I_current_all = [0.00];

%%
for iSet = 1 : length(I_current_all)
    J_min       = 0.21;  % (nA) the minimum recurrent connectivity strength
    J_max       = 0.3;  % (nA) the maximum recurrent connectivity strength
    G           = 0.48;
    
    % simulation setting
    Ttotal      = 15;  % Total duration of simulation (s)
    dt          = 0.01 / 1000;  % Simulation time step
    
    % stimulus input
    I_stim    = I_current_all(iSet); % (nA)
    
    % distractor input
    I_dist      = 0; % (nA)
    
    if exist('./data','dir' ) == 0
        mkdir('./data')
    end
    
    disp(['J_min=', num2str(J_min), '  J_max=', num2str(J_max), '  I_stim=', num2str(I_stim), '  I_dist=', num2str(I_dist)])
    
    if I_stim > 0 && I_dist == 0
        figure_folder_name = ['figure', '_', 'G', num2str(G), '_', 'Js', num2str(J_min), '-', num2str(J_max), '_', 'I', num2str(I_stim),  '_', 'Stim'];
        data_mat_name = ['./data/data', '_', 'G', num2str(G), '_', 'Js', num2str(J_min), '-', num2str(J_max),  '_', 'I', num2str(I_stim), '_', 'Stim', '.mat'];
    elseif I_stim == 0 && I_dist > 0
        figure_folder_name = ['figure', '_', 'G', num2str(G), '_', 'Js', num2str(J_min), '-', num2str(J_max), '_', 'I', num2str(I_dist),  '_', 'Dist'];
        data_mat_name = ['./data/data', '_', 'G', num2str(G), '_', 'Js', num2str(J_min), '-', num2str(J_max),'_', 'I', num2str(I_dist),  '_', 'Dist', '.mat'];
    elseif I_stim == 0 && I_dist == 0 
        figure_folder_name = ['figure', '_', 'G', num2str(G), '_', 'Js', num2str(J_min), '-', num2str(J_max), '_', 'I0'];
        data_mat_name = ['./data/data', '_', 'G', num2str(G), '_', 'Js', num2str(J_min), '-', num2str(J_max), '_', 'I0.mat'];
    elseif I_stim > 0 && I_dist > 0
        figure_folder_name = ['figure', '_', 'G', num2str(G), '_', 'Js', num2str(J_min), '-', num2str(J_max), '_', 'Istim', num2str(I_stim), '_', 'Idist', num2str(I_dist)];
        data_mat_name = ['./data/data', '_', 'G', num2str(G), '_', 'Js', num2str(J_min), '-', num2str(J_max),'_', 'Istim', num2str(I_stim), '_', 'Idist', num2str(I_dist), '.mat'];
    end 
    
    scatter_num = 10;
    r_end_all = zeros(scatter_num, area_num * 3);
    S_end_all = zeros(scatter_num, area_num * 3);
    init_all = zeros(scatter_num, area_num * 3);
    
    if exist(figure_folder_name, 'dir')
        rmdir(figure_folder_name, 's');
    end
    mkdir(figure_folder_name);
    
    for ii = 1 : scatter_num
        fprintf('processing: %d \n', ii)
        [init_value, r_end, S_end] = Sim_ode(ii, figure_folder_name, Ttotal, dt, J_min, J_max, G, 'I_stim', I_stim, 'I_dist', I_dist);
        init_all(ii, :) = init_value;
        r_end_all(ii, :) = r_end;
        S_end_all(ii, :) = S_end;
    end
    
    parsave(data_mat_name, r_end_all, S_end_all, init_all)
end
toc;

