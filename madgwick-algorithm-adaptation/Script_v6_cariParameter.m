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

%AHRS = MadgwickAHRS('SamplePeriod', 1/256, 'Beta', 0.02);
%AHRS = MahonyAHRS('SamplePeriod', 1/256, 'Kp', 0.9);

quaternion_sensor2earth = zeros(length(Accelerometer), 4);
betas = 0.001:0.0005:0.030;
gammas = 0:0.05:0.75;
% betas = 0.001;
% gammas = 0.1;
errors = zeros(length(betas), length(gammas));

for i = 1:length(betas)
    for j = 1:length(gammas)
        AHRS = MadgwickAHRS2('SamplePeriod', 1/256, 'Beta', betas(i), 'Gamma', gammas(j));
        for t = 1:length(Accelerometer)
            AHRS.Update(Gyroscope(t,:) * (pi/180), Accelerometer(t,:), Magnetometer(t,:));	% gyroscope units must be radians
            quaternion_sensor2earth(t, :) = quatconj(AHRS.Quaternion); % ini masih Earth relative to sensor
        end
        errors(i,j) = errorFunction(quaternion_sensor2earth, quaternionGT_sensor2earth);
    end
    disp(i);
end

