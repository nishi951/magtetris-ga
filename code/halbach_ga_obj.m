function [obj] = halbach_ga_obj( ...
    InnerR, array_spacing, R_diff_inout, ... % Optimization vars
    magnet_dim, Br, n_ring, half_FOV, FOV_1, FOV_2, FOV_3,... % Geometric parameters
    range_lambda, force_lambda, aniso_lambda, human_safe_force, use_force, verbose...           % Hyperparameters
    )
%HALBACH_GA_OBJ Summary of this function goes here
%   Detailed explanation goes here
[Halbach_MT, n_in_all, n_out_all] = get_mag_array( ...
    InnerR, array_spacing, R_diff_inout, ... % Optimization vars
    magnet_dim, Br, n_ring...           % Hyperparameters
);
n_per_group = 1;
obj = 0;

% Force penalty (very coarse)
min_ring = floor(n_ring/2); % Make sure this is >= 2
if use_force
    forces = halbach_force_calc(Halbach_MT, n_in_all, n_out_all, min_ring, 1);
    obj = obj + force_lambda * max(max(abs(forces(:))) - human_safe_force, 0);
end

for surface=['x', 'z']
    [By_mean, By_del] = field_mean_range(Halbach_MT, FOV_1, FOV_2, FOV_3, half_FOV, surface, n_per_group, verbose);

%     weight = Halbach_MT.Weight();
%     fprintf('Weight: %.1f kg\n', weight);
%     force = max(forces);
    obj = obj - aniso_lambda(surface) * By_mean;
    obj = obj + aniso_lambda(surface) * range_lambda * By_del;
end

end

