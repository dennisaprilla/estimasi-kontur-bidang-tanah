clear; clc; close all;

addpath('/home/dennis/Documents/MATLAB/dcm-algorithm-matlab-master/src');
addpath('/home/dennis/Documents/MATLAB/dcm-algorithm-matlab-master/quaternion_library');

load('quaternionGT_sensor2earth.mat');

%% Read data dan Bikin Noise

% Read IMU text
[Accelerometer, Gyroscope, Magnetometer] = readIMUtext('IMU.txt');

IMU = [Gyroscope, Accelerometer, Magnetometer];

% Bikin Noise multilevel
stationary_IMU = IMU(1:680,:);
stationary_IMU_normalized = stationary_IMU - mean(stationary_IMU, 1);
stationary_IMU_std = std(stationary_IMU_normalized, 0, 1);

noise_levels = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5];
% IMU_noised = zeros(size(IMU, 1), size(IMU, 2), length(noise_levels));
% for i = 1:length(noise_levels)
%     for axis = 1:size(IMU,2)
%         IMU_noised(:,axis,i) = IMU(:,axis) + normrnd(0, noise_levels(i)*stationary_IMU_std(axis), size(IMU(:,axis)));
%     end
% end


%% Mulai

% ceritanya time delta
sample_period = 1/256;

% buat nampung R matrix
quaternion_sensor2earth = zeros(length(Gyroscope), 4);

observations = 10;
DCMErrors = zeros(observations, length(noise_levels));

for observation = 1:observations
    
    fprintf('OBSERVATION %d ------------------ \n', observation);

    for noise_level = 1:length(noise_levels)
        
        IMU_noised = zeros(size(IMU));
        
        for axis = 1:size(IMU,2)
            IMU_noised(:,axis) = ...
                IMU(:,axis) + normrnd(0, noise_levels(noise_level)*stationary_IMU_std(axis), [size(IMU,1), 1]);
        end
        
        % init values of omega
        Omega_P = [0 0 0];
        Omega_I = [0 0 0];
        
        for i = 1:length(Gyroscope)
            if i==1
                % init the dcm matrix
                [ pitch, roll, yaw, dcm_init ] = reset_fusion( IMU_noised(i,4), ...
                                                               IMU_noised(i,5), ...
                                                               IMU_noised(i,6), ...
                                                               IMU_noised(i,7),  ...
                                                               IMU_noised(i,8),  ...
                                                               IMU_noised(i,9));
            else

                % update the dcm matrix with the previous one and the values
                [ pitch, roll, yaw, ...
                  dcm_init, ...
                  Omega_I,  ...
                  Omega_P ] = dcm_algorithm ( dcm_init, ...
                                              IMU_noised(i,4), ...
                                              IMU_noised(i,5), ...
                                              IMU_noised(i,6), ...
                                              IMU_noised(i,1)*(pi/180), ...
                                              IMU_noised(i,2)*(pi/180), ...
                                              IMU_noised(i,3)*(pi/180), ...
                                              IMU_noised(i,7), ...
                                              IMU_noised(i,8), ...
                                              IMU_noised(i,9), ...
                                              sample_period, ...
                                              pitch, roll, yaw, ...
                                              Omega_I, Omega_P );
            end

            quaternion_sensor2earth(i,:) = rotm2quat(dcm_init');
        end

        DCMErrors(observation, noise_level) = errorFunction(quaternion_sensor2earth, quaternionGT_sensor2earth);
        fprintf('noise level %d is done\n', noise_level);

    end
    
end

DCMErrors_meanObservation = mean(DCMErrors,1);

figure('Name', 'Errors');
grid on; hold on;
plot(noise_levels, DCMErrors_meanObservation, '-*r');

%save('DCMErrors.mat', 'DCMErrors_meanObservation', 'DCMErrors', 'noise_levels');