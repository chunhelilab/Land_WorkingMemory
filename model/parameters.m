clear;
clc;

%% area name and number
area_name = {'V1', 'V2', 'V4', 'DP', 'MT', '8m', '5', '8l', '2', 'TEO', 'F1', 'STPc', '7A', '46d', '10', '9/46v', '9/46d', 'F5', 'TEpd', 'PBr', '7m', 'LIP', 'F2', '7B', 'ProM', 'STPi', 'F7', '8B', 'STPr', '24c'};
area_num = 30;

%% parameters for simulation
tau_N       = 0.060;  % (s)
tau_G       = 0.005;  % (s)
tau_r       = 0.002; % (s)
gamma_E     = 1.282;
gamma_I     = 2;
I_0A        = 0.3294;  % (nA)
I_0B        = 0.3294;  % (nA)
I_0C        = 0.26;  % (nA)

%% transfer function for A and B population
a          = 135; % (Hz/nA)
b          = 54; % (Hz)
d          = 0.308; % (s)

%% transfer function for C population
g_I        = 4;
c_1        = 615; % (Hz/nA)
c_0        = 177; % (Hz)
r_0        = 5.5; % (Hz)

%% projection parameters for isolated area
J_s         = 0.3213; % (nA)
% J_c         = 0.0107; % (nA)
J_IE        = 0.15;  % (nA)
% J_EI        = -0.31;  % (nA)
% J_II        = -0.12;  % (nA)
% zeta = (tau_G * gamma_I * c_1) / (g_I - J_II * tau_G * gamma_I * c_1);
% J_0 = J_s + J_c + 2 * J_EI * J_IE * zeta;  % (nA)

%% Inter-area projection, SLN
SLN         = load("./SLN_matrix.mat").SLN;  % (to, from)
SLN_i = ones(area_num, area_num) - SLN;
frontal = [6 8 11 14 15 16 17 18 23 25 27 28]; % areas in frontal lobe
C_FEF = 0.4; % this means long-range inhibition to FEF tops at maxi
for i = 1 : 2 % 8m and 8l
    for j = 1 : length(frontal)
        x = frontal(i);
        y = frontal(j);
        if SLN_i(x,y) > C_FEF
            SLN_i(x,y) = C_FEF;
        end
    end
end

%% Inter-area projection, FLN
k_1         = 1.2;
k_2         = 0.3;
FLN         = load('./FLN_matrix.mat').FLN;  % (to, from)
FLN = k_1 .* FLN .^ (k_2);
norm = sum(FLN, 2);
for k = 1 : area_num
    FLN(k,:) = FLN(k,:) / norm(k, 1);
end

%% gradient of spine conuts
spine_count = load("./spine_count.mat").spine_count';
[~, spine_count_idx] = sort(spine_count);
area_name_ranked = area_name(1, spine_count_idx);

%% noise
tau_noise   = 0.002; 
sigma_A     = 0.0002; 
sigma_B     = 0.0002; 
sigma_C     = 0.0002;  

%% parameters for calculating transition paths
params = {};
params.N = 30; % The number of points in the minimum action path.
params.TMax = 10; % The time range of the minimum action path. Larger values can be more accurate, but can also lead to instabilities.
params.c = 1e10; % The remeshing parameter, c. Larger values give greater remeshing.
params.MaxIter = 1000; % The number of optimization steps to run between remeshing steps.
params.K = 3; % The number of total remeshing steps to do;
params.Remesh = 1; % Turn the remesh on or off. If off, the algorithm will be default perform K*MaxIter itereations.
params.q = 3; % The q parameter from the original paconstraints can be set of the path as well. Is a measure of path smoothness.
params.ub = []; % If desired, manuscript.
params.lb = []; % As above, lower bounds can be set on the path as well. This is also not done at all in our manuscript.
params.PhiInit =[]; % The default intial path between two states is a straight line. This can be modified here.
params.ObjGrad = 0; % params.ObjGrad = Parameters.JacobianPresent; We require a Jacobian in this implementation.
params.Dimension = 90; % The dimension of the system.

%%
save('parameters.mat')