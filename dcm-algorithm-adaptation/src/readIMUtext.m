function [ Accelerometer, Gyroscope, Magnetometer ] = readIMUtext( text )
%READIMUTEXT Summary of this function goes here
%   Detailed explanation goes here

    if(isa(text,'char'))
        % mempersiapkan parameter
        fileID = fopen(text,'r');
        formatSpec = '%f';
        sizeA = [9 Inf];
        
        % read file
        IMU = fscanf(fileID,formatSpec,sizeA)';
        
        % pecah imu
        Accelerometer = IMU(:, 1:3);
        Gyroscope     = IMU(:, 4:6);
        Magnetometer  = IMU(:, 7:9);
    end

end

