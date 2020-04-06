function [ pitch, roll, yaw, dcm_init ] = ...
         reset_fusion( accx, accy, accz, magx, magy, magz )
% reset_fusion  Computes the initial dcm matrix
%   Inputs:
%     -accx:  acceleration in the x axis [g]
%     -accy:  acceleration in the y axis [g]
%     -accz:  acceleration in the z axis [g]
%     -magx:  magnetic filed in the x axis [gauss]
%     -magy:  magnetic filed in the x axis [gauss]
%     -magz:  magnetic filed in the x axis [gauss]
%
%   Outputs:
%     -pitch:    pitch
%     -roll:     roll
%     -yaw:      yaw
%     -dcm_init: dcm matrix for the given values
%
%   See also:
%      -github.com/Razor-AHRS/razor-9dof-ahrs/tree/master/Arduino
%      -http://bth.diva-portal.org/smash/get/diva2:1127455/FULLTEXT02.pdf

% calculate pitch
pitch = -atan2(accx, sqrt(accy^2 + accz^2));

% calculate roll
temp1 = cross([accx accy accz], [1 0 0]);
temp2 = cross([1 0 0], temp1);
roll  = atan2(temp2(2), temp2(3));

% calculate yaw
mag_x = magx * cos(pitch) + magy * sin(roll) * sin(pitch) + magz * cos(roll) * sin(pitch);
mag_y = magy * cos(roll) - magz * sin(roll);
yaw   = atan2(-mag_y, mag_x);

% init rotation matrix
% this should have the same result
dcm_init = angle2dcm( pitch, roll, yaw, 'YXZ' );

% m00 = cos(pitch) * cos(yaw);
% m01 = cos(yaw)   * sin(roll) * sin(pitch) - cos(roll) * sin(yaw);
% m02 = sin(roll)  * sin(yaw) + cos(roll) * cos(yaw) * sin(pitch);
% 
% m10 = cos(pitch) * sin(yaw);
% m11 = cos(roll)  * cos(yaw) + sin(roll) * sin(pitch) * sin(yaw);
% m12 = cos(roll)  * sin(pitch) * sin(yaw) - cos(yaw) * sin(roll);
% 
% m20 = -sin(pitch);
% m21 = cos(pitch) * sin(roll);
% m22 = cos(roll)  * cos(pitch);
% 
% dcm_init = [m00 m01 m02; m10 m11 m12; m20 m21 m22];
end

