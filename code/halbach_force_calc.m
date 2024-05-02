function [F_targets] = halbach_force_calc(MT, n_in_all, n_out_all, min_ring, FOV_del_ratio)
n_ring = size(n_in_all, 2);
F_targets = zeros(3, n_ring-1);
out_idx = 1;
for target_layer_idx=min_ring:n_ring
    source_layer_idx = 1:n_ring-1;
    % Looping through target layer(s)
    num_target = sum(n_in_all(target_layer_idx)) + sum(n_out_all(target_layer_idx));
    target_idx_list = ones(num_target,2);
    target_end_idx = 0;
    for curr_layer_idx=1:length(target_layer_idx)
        curr_layer = target_layer_idx(curr_layer_idx);
        num_curr_layer = n_in_all(curr_layer) + n_out_all(curr_layer);
        start_idx = sum(n_in_all(1:(curr_layer-1))) + sum(n_out_all(1:(curr_layer-1)));
        target_start_idx = target_end_idx + 1;
        target_end_idx = target_end_idx + num_curr_layer;
        target_idx_list(target_start_idx:target_end_idx,2) = (start_idx+1):(start_idx+num_curr_layer);
    end
    
    % Looping through source layer(s)
    num_source = sum(n_in_all(source_layer_idx)) + sum(n_out_all(source_layer_idx));
    source_idx_list = ones(num_source,2);
    source_end_idx = 0;
    for curr_layer_idx=1:length(source_layer_idx)
        curr_layer = source_layer_idx(curr_layer_idx);
        num_curr_layer = n_in_all(curr_layer) + n_out_all(curr_layer);
        start_idx = sum(n_in_all(1:(curr_layer-1))) + sum(n_out_all(1:(curr_layer-1)));
        source_start_idx = source_end_idx + 1;
        source_end_idx = source_end_idx + num_curr_layer;
        source_idx_list(source_start_idx:source_end_idx,2) = (start_idx+1):(start_idx+num_curr_layer);
    end
    
    % Force calculation
    F_target_list = MT.ForceFree(target_idx_list,source_idx_list,FOV_del_ratio);
    F_targets(:, out_idx) = sum(F_target_list,1);
    out_idx = out_idx + 1;
end