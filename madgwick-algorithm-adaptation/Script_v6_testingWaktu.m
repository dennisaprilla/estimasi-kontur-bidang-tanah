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

%% Import and plot sensor data

%load('ExampleData.mat');
[Accelerometer, Gyroscope, Magnetometer] = readIMUtext('IMU.txt');


%% Process sensor data through algorithm

quaternion_sensor2earth = zeros(length(Accelerometer), 4);

beta = 0.001;
gamma = 0.1;
AHRS = MadgwickAHRS2('SamplePeriod', 1/256, 'Beta', beta, 'Gamma', gamma);

observations = 10;
madgwickDurationTimes = zeros(3, observations);
madgwickDurationTimes_mean = [];

%% Durasi hanya peritungan

for observation = 1:observations
    
    tic;
    for t = 1:length(Accelerometer)
        AHRS.Update(Gyroscope(t,:) * (pi/180), Accelerometer(t,:), Magnetometer(t,:));	% gyroscope units must be radians
    end
    madgwickDurationTimes(1,observation) = toc;
    
end

disp('Tahap 1 selesai');

%% Durasi dengan penyimpanan data

for observation = 1:observations
    
    tic;
    for t = 1:length(Accelerometer)
        AHRS.Update(Gyroscope(t,:) * (pi/180), Accelerometer(t,:), Magnetometer(t,:));	% gyroscope units must be radians
        quaternion_sensor2earth(t, :) = AHRS.Quaternion;
    end
    madgwickDurationTimes(2,observation) = toc;
    
end

disp('Tahap 2 selesai');

%% Durasi dengan konjugasi

for observation = 1:observations
    
    tic;
    for t = 1:length(Accelerometer)
        AHRS.Update(Gyroscope(t,:) * (pi/180), Accelerometer(t,:), Magnetometer(t,:));	% gyroscope units must be radians
        quaternion_sensor2earth(t, :) = quatconj(AHRS.Quaternion);
    end 
    madgwickDurationTimes(3,observation) = toc;
    
end

disp('Tahap 3 selesai');

%%

madgwickDurationTimes_mean = mean(madgwickDurationTimes,2)';

% save('madgwickDuration.mat', 'madgwickDurationTimes', 'madgwickDurationTimes_mean', 'quaternion_sensor2earth');