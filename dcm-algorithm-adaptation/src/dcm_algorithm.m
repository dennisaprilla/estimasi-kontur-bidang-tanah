function [ pitch_o, roll_o, yaw_o, dcm, Omega_I, Omega_P ] = ...
    dcm_algorithm( dcm_init, ...
                   accx, accy, accz, ...
                   gyrx, gyry, gyrz, ...
                   magx, magy, magz, ...
                   G_Dt, ...
                   pitch, roll, yaw, ...
                   Omega_I, Omega_P)
% dcm_algorithm  Computes and updates
%   Inputs:
%     -dcm_init: dcm matrix form last iteration
%     -accx:     acceleration in the x axis [g]
%     -accy:     acceleration in the y axis [g]
%     -accz:     acceleration in the z axis [g]
%     -gyrx:     angular speed in the x axis [deg/s]
%     -gyry:     angular speed in the y axis [deg/s]
%     -gyrz:     angular speed in the z axis [deg/s]
%     -magx:     magnetic filed in the x axis [gauss]
%     -magy:     magnetic filed in the x axis [gauss]
%     -magz:     magnetic filed in the x axis [gauss]
%     -G_Dt:     delta of time since last update [sec]
%     -pitch:    last pitch value [deg] 
%     -roll:     last roll value [deg]
%     -yaw:      last roll value [deg] (unused)
%     -Omega_I:  integral term for drift correction
%     -Omage_P:  proportional term for drift correction
%
%   Outputs:
%     -pitch_o: pitch
%     -roll_o:  roll
%     -yaw_o:   yaw
%     -dcm:     dcm matrix
%     -Omega_I: integral term for drift correction
%     -Omega_P: proportional term for drift correction
%
%   See also:
%      -github.com/Razor-AHRS/razor-9dof-ahrs/tree/master/Arduino
%      -http://bth.diva-portal.org/smash/get/diva2:1127455/FULLTEXT02.pdf

% tune this variables
% Kp_ROLLPITCH  = 0.02;
% Ki_ROLLPITCH  = 0.00002;
% Kp_YAW        = 1.2;
% Ki_YAW        = 0.00002;

% % best (untuk pengujian orientasi)
% Kp_ROLLPITCH  = 0.13;
% Ki_ROLLPITCH  = 0.0003;
% Kp_YAW        = 25;
% Ki_YAW        = 0.0083;

% % sensor sendiri (bagus untuk groundtruth 2, checkpoin 1)
% Kp_ROLLPITCH  = 0.1;
% Ki_ROLLPITCH  = 0;
% Kp_YAW        = 0.003;
% Ki_YAW        = 0.0000001;

% % sensor sendiri (bagus untuk groundtruth 2, checkpoin 2)
% Kp_ROLLPITCH  = 0.2;
% Ki_ROLLPITCH  = 0.00001;
% Kp_YAW        = 0.0001;
% Ki_YAW        = 0.0000001;

% % sensor sendiri (bagus untuk groundtruth 2, checkpoin 3)
% Kp_ROLLPITCH  = 1.8;
% Ki_ROLLPITCH  = 0.00001;
% Kp_YAW        = 0.0001;
% Ki_YAW        = 0.0000001;

% % sensor sendiri (bagus untuk groundtruth 2, checkpoin 4)
% Kp_ROLLPITCH  = 1.8;
% Ki_ROLLPITCH  = 0.00001;
% Kp_YAW        = 0.0000;
% Ki_YAW        = 0.0000001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % sensor sendiri (bagus untuk groundtruth 3, checkpoint 1)
% Kp_ROLLPITCH  = 1;
% Ki_ROLLPITCH  = 0.00001;
% Kp_YAW        = 0.0001;
% Ki_YAW        = 0.000001;

% % sensor sendiri (bagus untuk groundtruth 3, checkpoint 2)
% Kp_ROLLPITCH  = 0.1;
% Ki_ROLLPITCH  = 0.00001;
% Kp_YAW        = 0.00005;
% Ki_YAW        = 0.000001;

% % sensor sendiri (bagus untuk groundtruth 3, checkpoint 3)
% Kp_ROLLPITCH  = 0.1;
% Ki_ROLLPITCH  = -0.00002;
% Kp_YAW        = 0.00005;
% Ki_YAW        = 0.0000005;

% % sensor sendiri (bagus untuk groundtruth 3, checkpoint 4)
% Kp_ROLLPITCH  = 0.1;
% Ki_ROLLPITCH  = -0.00002;
% Kp_YAW        = 0.00005;
% Ki_YAW        = 0.0000005;

% % sensor sendiri (bagus untuk groundtruth 3 [bisa juga untuk 2], checkpoint 5)
% Kp_ROLLPITCH  = 0.5;
% Ki_ROLLPITCH  = -0.00001;
% Kp_YAW        = 0.0005;
% Ki_YAW        = 0.0000002;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % sensor sendiri (bagus untuk groundtruth 2, checkpoint 1)
% Kp_ROLLPITCH  = 0.02;
% Ki_ROLLPITCH  = 0.00001;
% Kp_YAW        = 0.01;
% Ki_YAW        = 0.00000002;

% sensor sendiri (bagus untuk groundtruth 2, checkpoint 5)
Kp_ROLLPITCH  = 0.1;
Ki_ROLLPITCH  = 0.00001;
Kp_YAW        = 0.00005;
Ki_YAW        = 0.0000005;

% % sensor sendiri (bagus untuk groundtruth 3, checkpoint 5)
% Kp_ROLLPITCH  = 0.02;
% Ki_ROLLPITCH  = -0.000004;
% Kp_YAW        = 0.0005;
% Ki_YAW        = 0.0000002;

% init the dcm matrix
Update_Matrix = [[]];

% Calculate yaw
mag_x = magx * cos(pitch) + ...
        magy * sin(roll) * sin(pitch) + ...
        magz * cos(roll) * sin(pitch);
mag_y = magy * cos(roll) - magz * sin(roll);
yaw = atan2(-mag_y, mag_x);

% Matrix update
Omega = [gyrx gyry gyrz] + Omega_I;
Omega_Vector = Omega + Omega_P;

% Use drift correction
Update_Matrix(1,1) = 0;
Update_Matrix(1,2) = -G_Dt*Omega_Vector(3); %-z
Update_Matrix(1,3) = G_Dt*Omega_Vector(2);  % y
Update_Matrix(2,1) = G_Dt*Omega_Vector(3);  % z
Update_Matrix(2,2) = 0;
Update_Matrix(2,3) = -G_Dt*Omega_Vector(1); %-x
Update_Matrix(3,1) = -G_Dt*Omega_Vector(2); %-y
Update_Matrix(3,2) = G_Dt*Omega_Vector(1);  % x
Update_Matrix(3,3) = 0;

% Matrix multiplication
Temp_Matrix = dcm_init * Update_Matrix;
dcm = Temp_Matrix + dcm_init;

% Normalize
temp = [[]];
error = -dot(dcm(1,:), dcm(2,:)) * .5;

% Multiply vector by scalar
% eq 19
temp(1,:) = dcm(2,:) * error;
temp(2,:) = dcm(1,:) * error;

temp(1,:) = temp(1,:) + dcm(1,:);
temp(2,:) = temp(2,:) + dcm(2,:);

% eq 20
temp(3,:) = cross(temp(1,:), temp(2,:));

% % eq 21
% renorm = .5 * ( 3 - dot(temp(1,:), temp(1,:)));
% dcm(1,:) = temp(1,:) * renorm;
% 
% renorm = .5 * ( 3 - dot(temp(2,:), temp(2,:)));
% dcm(2,:) = temp(2,:) * renorm;
% 
% renorm = .5 * ( 3 - dot(temp(3,:), temp(3,:)));
% dcm(3,:) = temp(3,:) * renorm;

dcm(1,:) = temp(1,:) / norm(temp(1,:));
dcm(2,:) = temp(2,:) / norm(temp(2,:));
dcm(3,:) = temp(3,:) / norm(temp(3,:));


% drift correction
accel_magnitude = sqrt(accx^2 + accy^2 + accz^2);

% roll and pitch
accel_weight = 1 - 2*abs( 1 - accel_magnitude );

% limit to range 0-1
accel_weight(accel_weight > 1) = 1;
accel_weight(accel_weight < 0) = 0;

errorRollPitch = cross ([accx accy accz], dcm(3,:));
Omega_P = errorRollPitch * Kp_ROLLPITCH*accel_weight;

Scaled_Omega_I = errorRollPitch * Ki_ROLLPITCH*accel_weight;
Omega_I = Omega_I + Scaled_Omega_I;

% yaw
mag_heading_x = cos(yaw);
mag_heading_y = sin(yaw);

errorCourse=(dcm(1,1)*mag_heading_y) - (dcm(2,1)*mag_heading_x);
errorYaw = dcm(3,:) * errorCourse;

Scaled_Omega_P = errorYaw * Kp_YAW;
Omega_P = Omega_P + Scaled_Omega_P;

Scaled_Omega_I = errorYaw * Ki_YAW;
Omega_I = Omega_I + Scaled_Omega_I;

% convert to euler angles
pitch_o = (-asin(dcm(3,1)));
roll_o  = (atan2(dcm(3,2),dcm(3,3)));
yaw_o   = (atan2(dcm(2,1),dcm(1,1)));

end
