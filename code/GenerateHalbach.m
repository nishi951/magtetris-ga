function [loc_all,angle_all,magnet_dim_all,Br_all] = GenerateHalbach(n_bar,R_in,magnet_dim,Br,angle_array_rot,z_center)
% n_bar: Number of magnets in this ring
% R_in: Radius of the ring
% magnet_dim: Dimensions of the magnet
% Br: Remanance of each magnet
% angle_array_rot: Global rotation of the ring
% z_center: z_offset of the ring (defaults to 0)

% ===== Array definition =====
angle_all_orig = zeros(n_bar,3);
magnet_dim_all = repmat(magnet_dim,n_bar,1);
Br_all = ones(1,n_bar)*Br;

% Angle
angle_del = -180/(n_bar/4);
angle_yaw_curr = (0:1:(n_bar-1))*angle_del;
angle_yaw_curr = angle_yaw_curr';
angle_all_orig(:,1) = angle_yaw_curr;

% Location
R_real_curr = R_in*ones(1,n_bar) + magnet_dim(2)/2;
loc_all_xy = zeros(n_bar,2);
for i=1:n_bar
    loc_curr_start_xy = [0,R_real_curr(i)]';
    angle_curr = -360/n_bar*(i-1);
    % Rotation about z-axis
    rot_mat = [cosd(angle_curr),-sind(angle_curr);...
               sind(angle_curr),cosd(angle_curr)];
    loc_all_xy(i,:) = (rot_mat*loc_curr_start_xy)';
end
if nargin == 6
    loc_all_z = ones(n_bar,1)*z_center;
else
    loc_all_z = zeros(n_bar,1);
end
loc_all_orig = [loc_all_xy,loc_all_z];

% Rotate the whole array
[loc_all,angle_all] = RotateArray(loc_all_orig,angle_all_orig,angle_array_rot,1);

end

