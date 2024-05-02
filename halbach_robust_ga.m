%% Define Variables
clear all;
close all;
clc;

verbose = false;
addpath ./code

% Constants
Br_N48 = 1.38;  % Remanance (NdFeB, N48)
% Br_N35 = ...  % TODO: Add other remanences...
FOV_del_ratio = 20; % Determines accuracy of force calculation

% Hyperparameters (change)
Br = Br_N48;
magnet_dim = [12, 12, 12]; % NdFeB magnet dimensions
half_FOV = 10; % 1/2 FOV, [mm]
human_safe_force = 144; % Target maximum force between planes, [N]
n_ring = 6;  % Number of z planes of rings
min_R_diff = 25; % Minimum spacing between inner and outer ring
min_spacing = 24; % 2x width of magnet
max_spacing = 32; % Reduce search space [mm]
range_lambda = 1e8; % Loss weight for field inhomogeneity
aniso_lambda = dictionary('x', 5., 'z', 1.); % Relative weights for homogeneity in different planes
force_lambda = 100; % Loss weight for applied force above human-safe force
n_diff_inout = 4; % Difference between inner and outer layers
use_force = false;
% max_per_ring = 30; % Maximum number of magnets per ring

% Parameters for robust genetic algorithm
n_samples = 3;
Br_std = 1e-2;
robust_mode = 'worst';

% Derived parameters
dx_opt = 5;
FOV_1 = -half_FOV:dx_opt:half_FOV;
FOV_2 = -half_FOV:dx_opt:half_FOV;
FOV_3 = 0; % Isocenter
% R_diff = min_R_diff; % Fix to be optimal

% Optimization variables (only need to optimize half the rings and flip)
min_R = (3*half_FOV + min_R_diff/2) .* ones(1, n_ring/2);
min_N = zeros(1, n_ring/2);
% InnerR = optimvar("InnerR", n_ring/2, "Type", "continuous", "LowerBound", min_R);   % Inner radius
% InnerN = optimvar("InnerN", n_ring/2, "Type", "integer", "LowerBound", min_N);
% OuterR = optimvar("OuterR", n_ring/2, "Type", "continuous", "LowerBound", min_R);   % Outer radius
% OuterN = optimvar("OuterN", n_ring/2, "Type", "integer", "LowerBound", min_N);
% R_diff = optimvar("R_diff", "Type", "continuous", "LowerBound", min_R_diff);  % Difference in radius from inner to outer, [mm]
% array_spacing = optimvar("array_spacing", n_ring/2, "Type", "continuous", "LowerBound", zeros(1, n_ring/2), "UpperBound", max_spacing*ones(1, n_ring/2)); % Layer spacing

% Test fitness function
% InnerR = [10.; 20.; 30.];
% InnerN = [13; 13; 13];
% array_spacing = [14.; 14.; 14;];
% n_ring = 6;
% halbach_ga_obj( ...
%     InnerR, array_spacing, ...
%     min_R_diff, n_diff_inout, magnet_dim, Br, n_ring, half_FOV, FOV_1, FOV_2, FOV_3, lambda...
% );

% Optimization objective
% Mean field
% Variable structure:
% - InnerR - continuous  [n_ring/2]
% - Array Spacing - continuous [n_ring/2]
% - R_diff  [n_ring / 2]
fitnessfcn = @(vars)( ...
    robust_halbach_ga_obj( ...
    vars(1:n_ring/2), ...         % InnerR
    vars(n_ring/2+1:n_ring),  ...    % array_spacing
    vars(n_ring+1:end), ...       % R_diff
    magnet_dim, Br, n_ring, half_FOV, FOV_1, FOV_2, FOV_3, ...
    robust_mode, Br_std, n_samples, ...
    range_lambda, force_lambda, aniso_lambda, human_safe_force, use_force, verbose...
));


% Constraints
InnerR_lb = (half_FOV + min_R_diff/2) * ones(n_ring/2, 1);
InnerR_ub = Inf(n_ring/2, 1);
% InnerR_ub = InnerR_lb + 1.;
array_spacing_lb = min_spacing * ones(n_ring/2, 1);
array_spacing_ub = max_spacing * ones(n_ring/2, 1);
R_diff_lb = min_R_diff * ones(n_ring/2, 1);
R_diff_ub = Inf(n_ring/2, 1);
% R_diff_ub = (min_R_diff + 0.1) * ones(n_ring/2, 1);
lb = [InnerR_lb; array_spacing_lb; R_diff_lb];
ub = [InnerR_ub; array_spacing_ub; R_diff_ub];


%% Solve (first stage, no force constraint)
opts = optimoptions(@ga, ...
                    'PopulationSize', 3000, ...
                    'MaxGenerations', 30, ...
                    'EliteCount', 60, ...
                    'FunctionTolerance', 1e-8, ...
                    'PlotFcn', @gaplotbestf, ...
                    'Display', 'diagnose');

rng default % For reproducibility
[sol,fval,exitflag,output,population,scores] = ga(fitnessfcn, 3*n_ring/2,[], [], [], [], lb, ub, [], [], opts);
save("population_noforce.mat", "population")

%% Solve (second stage, force constraint)
% use_force = true;
% opts = optimoptions(@ga, ...
%                     'PopulationSize', 40, ...
%                     'MaxGenerations', 4, ...
%                     'EliteCount', 4, ...
%                     'FunctionTolerance', 1e-8, ...
%                     'PlotFcn', @gaplotbestf);
% opts.InitialPopulation = population;
% 
% fitnessfcn = @(InnerR_array_spacing)( ...
%     halbach_ga_obj( ...
%     InnerR_array_spacing(1:n_ring/2), InnerR_array_spacing(n_ring/2+1:end), ...
%     min_R_diff, n_diff_inout, magnet_dim, Br, n_ring, half_FOV, FOV_1, FOV_2, FOV_3, ...
%     range_lambda, force_lambda, aniso_lambda, human_safe_force, use_force...
% ));
% [sol,fval,exitflag,output,population,scores] = ga(fitnessfcn, n_ring,[], [], [], [], lb, ub, [], [], opts);


%% Analyze
InnerR_sol = sol(1:n_ring/2);
array_spacing_sol = sol(n_ring/2+1:n_ring);
R_diff_sol = sol(n_ring+1:end);
% InnerR_sol = population(1, 1:n_ring/2);
% array_spacing_sol = population(1, n_ring/2+1:n_ring);
% R_diff_sol
[MT, n_in_all, n_out_all] = get_mag_array(InnerR_sol, array_spacing_sol, R_diff_sol, magnet_dim, Br, n_ring);

%% Display
dx_val = 1.;
FOV_1 = -half_FOV:dx_val:half_FOV;
FOV_2 = -half_FOV:dx_val:half_FOV;
FOV_3 = 0; % Isocenter
surface = 'x';
[Bx, By, Bz] = MT.Field2D(FOV_1, FOV_2, FOV_3, surface, 1);
% Only keep the circular region
for idx_1=1:length(FOV_1)
    for idx_2=1:length(FOV_2)
        r = sqrt(FOV_1(idx_1)^2 + FOV_2(idx_2)^2);
        if r > half_FOV
            Bx(idx_2,idx_1) = NaN;
            By(idx_2,idx_1) = NaN;
            Bz(idx_2,idx_1) = NaN;
        end
    end
end
Bx = Bx * 1e3;
By = By * 1e3;
Bz = Bz * 1e3;

DrawMagnetization(MT.cub_loc, MT.cub_angle, MT.cub_dim, FOV_1, FOV_2, FOV_3, surface, By);

%% Compute final Mean field and inhomogeneity
% Derived parameters

[By_mean_final, By_del_final] = field_mean_range(MT, FOV_1, FOV_2, FOV_3, half_FOV, surface, 1, verbose);
fprintf('Final Bandwidth: %.3f%%\n', By_del_final/By_mean_final*100);
fprintf('Final Mean: %.1f mT\n', By_mean_final);
%% Compute maximum force
forces = halbach_force_calc(MT, n_in_all, n_out_all, 1, 20);
fprintf('Maximum force: %.3f \n', max(abs(forces(:))));

%% Save
save("MT_final.mat", "MT", "n_in_all", "n_out_all", "sol");
