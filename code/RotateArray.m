function [loc_final,angle_final] = RotateArray(loc_orig,angle_orig,angle_rot,rot_mode,rot_center)
% Rotate the whole array by angle_rot(yaw,pitch,row)
% INPUT:
%       loc_orig - the original center locations of the magnets
%       angle_orig - the original orientation angles of the magnets
%       angle_rot - the angle to rotate the whole array
%       rot_mode - 1: roll-pitch-yaw, 2: yaw-pitch-roll
% OUTPUT:
%       loc_final - the final center locations of the magnets
%       angle_final - the final orientation angles of the magnets

angle_final = angle_orig;
angle_final(:,1) = angle_final(:,1) + angle_rot(1);
angle_final(:,2) = angle_final(:,2) + angle_rot(2);
angle_final(:,3) = angle_final(:,3) + angle_rot(3);

angle_yaw = angle_rot(1);
angle_pitch = angle_rot(2);
angle_roll = angle_rot(3);
rot_roll = [1,0,0;...
            0,cosd(angle_roll),-sind(angle_roll);...
            0,sind(angle_roll),cosd(angle_roll)];
rot_pitch = [cosd(angle_pitch),0,sind(angle_pitch);...
             0,1,0;...
             -sind(angle_pitch),0,cosd(angle_pitch)];
rot_yaw = [cosd(angle_yaw),-sind(angle_yaw),0;...
           sind(angle_yaw),cosd(angle_yaw),0;...
           0,0,1];
% Choose different rotation modes
if rot_mode == 1
    rot_mat = (rot_roll*rot_pitch*rot_yaw)';
elseif rot_mode == 2
    rot_mat = (rot_yaw*rot_pitch*rot_roll)';
end
loc_final = loc_orig*rot_mat;
end

