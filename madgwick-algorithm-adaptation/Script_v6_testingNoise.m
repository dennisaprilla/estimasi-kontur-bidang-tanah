% ExampleScript.m
%
% This script demonstrates use of the MadgwickAHRS and MahonyAHRS algorithm
% classes with example data. ExampleData.mat contains calibrated gyroscope,
% accelerometer and magnetometer data logged from an AHRS device (x-IMU)
% while it was sequentially rotated from 0 degrees, to +90 degree and then
% to -90 degrees around the X, Y and Z axis.  The script first plots the
% example sensor data, then processes the data through the algorithm and
% plots the output as Euler angles.
%
% Note that the Euler angle plot shows erratic behaviour in phi and psi
% when theta approaches �90 degrees. This due to a singularity in the Euler
% angle sequence known as 'Gimbal lock'.  This issue does not exist for a
% quaternion or rotation matrix representation.
%
% Date          Author          Notes
% 28/09/2011    SOH Madgwick    Initial release
% 13/04/2012    SOH Madgwick    deg2rad function no longer used
% 06/11/2012    Seb Madgwick    radian to degrees calculation corrected

%% Start of script

close all;                          % close all figures
clear;                              % clear all variables
clc;                                % clear the command terminal

load('quaternionGT_sensor2earth.mat');
addpath('quaternion_library'); 

%% Import and plot sensor data

%load('ExampleData.mat');
[Accelerometer, Gyroscope, Magnetometer] = readIMUtext('IMU.txt');

IMU = [Gyroscope, Accelerometer, Magnetometer];

%% Bikin Noise multilevel

stationary_IMU = IMU(1:680,:);
stationary_IMU_normalized = stationary_IMU - mean(stationary_IMU, 1);
stationary_IMU_std = std(stationary_IMU_normalized, 0, 1);

noise_levels = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5];
% IMU_noised = zeros(size(IMU, 1), size(IMU, 2), length(noise_levels));
% for noise_level = 1:length(noise_levels)
%     for j = 1:size(IMU,2)
%         IMU_noised(:,j,noise_level) = IMU(:,j) + normrnd(0, noise_levels(noise_level)*stationary_IMU_std(j), size(IMU(:,j)));
%     end
% end


%% Process sensor data through algorithm

beta = 0.001;
gamma = 0.1;

quaternion_sensor2earth = zeros(length(Gyroscope), 4);

observations = 10;
madgwickErrors = zeros(observations, length(noise_levels));

for observation = 1:observations
    
    fprintf('OBSERVATION %d ------------------ \n', observation);
    
    AHRS = MadgwickAHRS2('SamplePeriod', 1/256, 'Beta', beta, 'Gamma', gamma);
    
    for noise_level=1:length(noise_levels)
        
        IMU_noised = zeros(size(IMU));
        
        for axis = 1:size(IMU,2)
            IMU_noised(:,axis) = ...
                IMU(:,axis) + normrnd(0, noise_levels(noise_level)*stationary_IMU_std(axis), [size(IMU,1), 1]);
        end
        
        for t = 1:length(Gyroscope)

            AHRS.Update(IMU_noised(t, 1:3) * (pi/180), IMU_noised(t, 4:6), IMU_noised(t, 7:9));
            quaternion_sensor2earth(t, :) = quatconj(AHRS.Quaternion); % ini masih Earth relative to sensor

        end

        madgwickErrors(observation, noise_level) = errorFunction(quaternion_sensor2earth, quaternionGT_sensor2earth);
        fprintf('noise level %d is done\n', noise_level);
        
    end
end

madgwickErrors_meanObservation = mean(madgwickErrors,1);

figure('Name', 'Errors');
grid on; hold on;
plot(noise_levels, madgwickErrors_meanObservation, '-*r');

%save('madgwickErrors.mat', 'madgwickErrors_meanObservation', 'madgwickErrors', 'noise_levels');