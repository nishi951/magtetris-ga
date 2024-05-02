%% Analyze
clear all; close all;
verbose = true;
load("population_noforce.mat");
Br = 1.38;
magnet_dim = [12, 12, 12]; % NdFeB magnet dimensions
n_ring = 6;

sol = population(1500, :);
InnerR_sol = sol(1:n_ring/2);
array_spacing_sol = sol(n_ring/2+1:n_ring);
R_diff_sol = sol(n_ring+1:end);
% InnerR_sol = population(1, 1:n_ring/2);
% array_spacing_sol = population(1, n_ring/2+1:n_ring);
% R_diff_sol
[MT, n_in_all, n_out_all] = get_mag_array(InnerR_sol, array_spacing_sol, R_diff_sol, magnet_dim, Br, n_ring);

%% Display
dx_val = 1.;
half_FOV = 10; %[mm]
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
save("MT_population.mat", "MT", "n_in_all", "n_out_all", "sol");
