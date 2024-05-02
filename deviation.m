%% Compute final Mean field and inhomogeneity
clear all; close all;

load('MT_robust.mat');
dx_val = 1.;
half_FOV = 10; % [mm]
verbose = false;
FOV_1 = -half_FOV:dx_val:half_FOV;
FOV_2 = -half_FOV:dx_val:half_FOV;
FOV_3 = 0; % Isocenter
surface = 'x';
[Bx, By, Bz] = MT.Field2D(FOV_1, FOV_2, FOV_3, surface, 1);
% DrawMagnetization(MT.cub_loc, MT.cub_angle, MT.cub_dim, FOV_1, FOV_2, FOV_3, surface, By);

[By_mean_final, By_del_final] = field_mean_range(MT, FOV_1, FOV_2, FOV_3, half_FOV, surface, 1, verbose);
fprintf('Final Bandwidth: %.3f%%\n', By_del_final/By_mean_final*100);
fprintf('Final Mean: %.1f mT\n', By_mean_final);

% Perturbation
n_trials = 100;
average_B0_robust = [];
average_B0_del_robust = [];
std_devs = logspace(-8, -2, 20);
for std_dev = std_devs
    means = [];
    bandwidths = [];
    for trial = 1:n_trials
        disp(trial);
        nominalBr = MT.cub_Br(1);
        deltaBr = std_dev * randn(size(MT.cub_Br));
        MT.cub_Br = MT.cub_Br + deltaBr;
        [By_mean_final, By_del_final] = field_mean_range(MT, FOV_1, FOV_2, FOV_3, half_FOV, surface, 1, verbose);
        means = cat(1, means, By_mean_final);
        bandwidths = cat(1, bandwidths, By_del_final);
    %     fprintf('Br_nominal: %.3\n', nominalBr);
    %     fprintf('Br_std_dev: %.3\n', std_dev);
    %     fprintf('Perturbed Final Bandwidth: %.3f%%\n', By_del_final/By_mean_final*100);
    %     fprintf('Perturbed Final Mean: %.1f mT\n', By_mean_final);
    end
    average_B0_robust = cat(1, average_B0_robust, mean(means));
    average_B0_del_robust = cat(1, average_B0_del_robust, mean(bandwidths));
end
f1 = figure;
loglog(std_devs, average_B0_robust, "DisplayName", "Robust");
title(sprintf('Average B0 (n=%d)', n_trials));
xlabel('standard dev of Br');
ylabel('field [mT]');
legend;
hold on;
f2 = figure;
loglog(std_devs, average_B0_del_robust, "DisplayName", "Robust");
title(sprintf('Average B0 delta (n=%d)', n_trials));
xlabel('standard dev of Br');
ylabel('inhomogeneity [%]');
legend;
hold on;
%%
load('MT_nonrobust.mat');
dx_val = 1.;
half_FOV = 10; % [mm]
verbose = false;
FOV_1 = -half_FOV:dx_val:half_FOV;
FOV_2 = -half_FOV:dx_val:half_FOV;
FOV_3 = 0; % Isocenter
surface = 'x';
[Bx, By, Bz] = MT.Field2D(FOV_1, FOV_2, FOV_3, surface, 1);
% DrawMagnetization(MT.cub_loc, MT.cub_angle, MT.cub_dim, FOV_1, FOV_2, FOV_3, surface, By);

[By_mean_final, By_del_final] = field_mean_range(MT, FOV_1, FOV_2, FOV_3, half_FOV, surface, 1, verbose);
fprintf('Final Bandwidth: %.3f%%\n', By_del_final/By_mean_final*100);
fprintf('Final Mean: %.1f mT\n', By_mean_final);

% Perturbation
n_trials = 100;
average_B0 = [];
average_B0_del = [];
std_devs = logspace(-8, -2, 20);
for std_dev = std_devs
    means = [];
    bandwidths = [];
    for trial = 1:n_trials
        disp(trial);
        nominalBr = MT.cub_Br(1);
        deltaBr = std_dev * randn(size(MT.cub_Br));
        MT.cub_Br = MT.cub_Br + deltaBr;
        [By_mean_final, By_del_final] = field_mean_range(MT, FOV_1, FOV_2, FOV_3, half_FOV, surface, 1, verbose);
        means = cat(1, means, By_mean_final);
        bandwidths = cat(1, bandwidths, By_del_final);
    %     fprintf('Br_nominal: %.3\n', nominalBr);
    %     fprintf('Br_std_dev: %.3\n', std_dev);
    %     fprintf('Perturbed Final Bandwidth: %.3f%%\n', By_del_final/By_mean_final*100);
    %     fprintf('Perturbed Final Mean: %.1f mT\n', By_mean_final);
    end
    average_B0 = cat(1, average_B0, mean(means));
    average_B0_del = cat(1, average_B0_del, mean(bandwidths));
end
figure(f1);
loglog(std_devs, average_B0, "DisplayName", "Non-robust");
% title(sprintf('Non-Robust Average B0 (n=%d)', n_trials));
xlabel('standard dev of Br');
ylabel('field [mT]');
legend;
figure(f2);
loglog(std_devs, average_B0_del, "DisplayName", "Non-robust");
% title(sprintf('Non-Robust Average B0 delta (n=%d)', n_trials));
xlabel('standard dev of Br');
ylabel('inhomogeneity [%]');

legend;

%%
load('MT_population.mat');
dx_val = 1.;
half_FOV = 10; % [mm]
verbose = false;
FOV_1 = -half_FOV:dx_val:half_FOV;
FOV_2 = -half_FOV:dx_val:half_FOV;
FOV_3 = 0; % Isocenter
surface = 'x';
[Bx, By, Bz] = MT.Field2D(FOV_1, FOV_2, FOV_3, surface, 1);
% DrawMagnetization(MT.cub_loc, MT.cub_angle, MT.cub_dim, FOV_1, FOV_2, FOV_3, surface, By);

[By_mean_final, By_del_final] = field_mean_range(MT, FOV_1, FOV_2, FOV_3, half_FOV, surface, 1, verbose);
fprintf('Final Bandwidth: %.3f%%\n', By_del_final/By_mean_final*100);
fprintf('Final Mean: %.1f mT\n', By_mean_final);

% Perturbation
n_trials = 100;
average_B0_pop = [];
average_B0_del_pop = [];
std_devs = logspace(-8, -2, 20);
for std_dev = std_devs
    means = [];
    bandwidths = [];
    for trial = 1:n_trials
        disp(trial);
        nominalBr = MT.cub_Br(1);
        deltaBr = std_dev * randn(size(MT.cub_Br));
        MT.cub_Br = MT.cub_Br + deltaBr;
        [By_mean_final, By_del_final] = field_mean_range(MT, FOV_1, FOV_2, FOV_3, half_FOV, surface, 1, verbose);
        means = cat(1, means, By_mean_final);
        bandwidths = cat(1, bandwidths, By_del_final);
    %     fprintf('Br_nominal: %.3\n', nominalBr);
    %     fprintf('Br_std_dev: %.3\n', std_dev);
    %     fprintf('Perturbed Final Bandwidth: %.3f%%\n', By_del_final/By_mean_final*100);
    %     fprintf('Perturbed Final Mean: %.1f mT\n', By_mean_final);
    end
    average_B0_pop = cat(1, average_B0_pop, mean(means));
    average_B0_del_pop = cat(1, average_B0_del_pop, mean(bandwidths));
end
figure(f1);
loglog(std_devs, average_B0_pop, "DisplayName", "Population");
% title(sprintf('Non-Robust Average B0 (n=%d)', n_trials));
xlabel('standard dev of Br');
ylabel('field [mT]');
legend;
figure(f2);
loglog(std_devs, average_B0_del_pop, "DisplayName", "Population");
% title(sprintf('Non-Robust Average B0 delta (n=%d)', n_trials));
xlabel('standard dev of Br');
ylabel('inhomogeneity [%]');

legend;
