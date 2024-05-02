function [obj] = robust_halbach_ga_obj( ...
    InnerR, array_spacing, R_diff_inout, ... % Optimization vars
    magnet_dim, Br, n_ring, half_FOV, FOV_1, FOV_2, FOV_3,... % Geometric parameters
    robust_mode, Br_std, n_samples, ...                       % Robust optimization params
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

field_objs = [];
nominal_Br = Halbach_MT.cub_Br;
for trial=1:n_samples
    field_obj = 0;
    Halbach_MT.cub_Br = nominal_Br + Br_std * randn(size(nominal_Br));
    for surface=['x', 'z']
        [By_mean, By_del] = field_mean_range(Halbach_MT, FOV_1, FOV_2, FOV_3, half_FOV, surface, n_per_group, verbose);
    
    %     weight = Halbach_MT.Weight();
    %     fprintf('Weight: %.1f kg\n', weight);
    %     force = max(forces);
%         field_obj = field_obj - aniso_lambda(surface) * By_mean;
        field_obj = field_obj - By_mean;
        field_obj = field_obj + aniso_lambda(surface) * range_lambda * By_del;
    end
    field_objs = cat(1, field_objs, field_obj);
end
if robust_mode == "mean"
    obj = obj + mean(field_objs);
elseif robust_mode == "worst"
    obj = obj + max(field_objs);
else % Just use the first one
    obj = obj + field_objs(1);
end

end

