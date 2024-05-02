function [MT, n_in_all, n_out_all] = get_mag_array( ...
    InnerR, array_spacing, R_diff_inout,... % Optimization vars
    magnet_dim, Br, n_ring ...           % Hyperparameters
 )
% Generate 2 layers of Halbach array for each ring (the other half is symmetric)
% Compute the number of magnets in each ring
angle_array_rot = [0, 0, 0];
% n_ring = 2*size(InnerR, 1);

% n_ring = 10;
% array_spacing = [1; 2; 3; 4; 5];
% Set up array spacing
if ~mod(n_ring, 2)
    % Even # of rings - Should have n_ring/2 spacings
    spacing = [array_spacing, array_spacing(end-1:-1:1)];
    z_centers = [0, cumsum(spacing)];
    z_centers = z_centers - z_centers(n_ring/2) - array_spacing(end) / 2;
else
    % Odd # of rings Should have n_ring/2 + 1 spacings
    spacing = [array_spacing, array_spacing(end:-1:1)];
    z_centers = [0, cumsum(spacing)];
    z_centers = z_centers - z_centers((n_ring+1)/2);
end
%%

R_in_all = zeros(1,n_ring);
R_out_all = zeros(1,n_ring);
n_in_all = zeros(1,n_ring);
n_out_all = zeros(1,n_ring);
for ring_idx=1:(n_ring/2)
    R_in_curr = InnerR(ring_idx);
    R_out_curr = R_in_curr + R_diff_inout(ring_idx);
%     n_in_curr = ceil(R_in_curr/3);
    n_in_curr = ceil(R_in_curr*2/5);
%     n_out_curr = n_in_curr + n_diff_inout;
%     n_out_curr = ceil(R_out_curr/3);
    n_out_curr = ceil(R_out_curr*2/5);
    % Get the symmetric index
    idx_sym = n_ring + 1 - ring_idx;
    R_in_all(ring_idx) = R_in_curr;
    R_out_all(ring_idx) = R_out_curr;
    n_in_all(ring_idx) = n_in_curr;
    n_out_all(ring_idx) = n_out_curr;
    R_in_all(idx_sym) = R_in_curr;
    R_out_all(idx_sym) = R_out_curr;
    n_in_all(idx_sym) = n_in_curr;
    n_out_all(idx_sym) = n_out_curr;
end
% Total number of magnets
n_total = sum(n_in_all) + sum(n_out_all);
loc_all = zeros(n_total,3);
angle_all = zeros(n_total,3);
Br_all = zeros(1,n_total);
magnet_dim_all = zeros(n_total,3);

% Generate the Halbach arrays
end_idx = 0;
for ring_idx=1:n_ring
    % z-location of the center of the current ring
%     z_center = (ring_idx-1 - (n_ring-1)/2)*(array_spacing(ring_idx) + magnet_dim(3));
    z_center = z_centers(ring_idx);
    % Inner ring
    start_idx = end_idx + 1;
    end_idx = end_idx + n_in_all(ring_idx);
    n_bar = n_in_all(ring_idx);
    R_inner = R_in_all(ring_idx);
    [loc_curr,angle_curr,magnet_dim_curr,Br_curr] = GenerateHalbach(n_bar,R_inner,magnet_dim,Br,angle_array_rot,z_center);
    loc_all(start_idx:end_idx,:) = loc_curr;
    angle_all(start_idx:end_idx,:) = angle_curr;
    magnet_dim_all(start_idx:end_idx,:) = magnet_dim_curr;
    Br_all(start_idx:end_idx) = Br_curr;
    % Outer ring
    start_idx = end_idx + 1;
    end_idx = end_idx + n_out_all(ring_idx);
    n_bar = n_out_all(ring_idx);
    R_inner = R_out_all(ring_idx);
    [loc_curr,angle_curr,magnet_dim_curr,Br_curr] = GenerateHalbach(n_bar,R_inner,magnet_dim,Br,angle_array_rot,z_center);
    loc_all(start_idx:end_idx,:) = loc_curr;
    angle_all(start_idx:end_idx,:) = angle_curr;
    magnet_dim_all(start_idx:end_idx,:) = magnet_dim_curr;
    Br_all(start_idx:end_idx) = Br_curr;
end

MT = MagTetris();
MT = MT.AssignCuboid(loc_all,angle_all,Br_all,magnet_dim_all);
end
