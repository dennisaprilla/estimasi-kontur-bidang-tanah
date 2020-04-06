% this file provides a quick test for verifying that the dcm matlab
% implementation gives the same result as the one implemented in c
% that can be found in [1].

% if you want to reproduce the same test, you will need to have the
% hardware setup described in [1] connected to the computer via a serial
% usb connexion, and senting data as a string with the following encoding:
%
%     [time, accx, accy, accz, gyrx, gyry, gyrz, magx, magy, magz ,...
%  ...expected_pitch, expected_roll, expecteD_yaw]

% note that the expected values are the ones calculated with the dcm
% implementation in c by the arduino, and are the values that we will
% use as a referenfe to verify ours,

% [1] github.com/Razor-AHRS/razor-9dof-ahrs/tree/master/Arduino

clear all; clc; close all;

try
    fclose(instrfindall);
catch
end

% modify the port number
% windows will need to change this
serialPort = serial('/dev/cu.usbmodem1411');
fopen(serialPort);

% init the arrays
Accelerometer = [[]];
Gyroscope     = [[]];
Magnetometer  = [[]];
Euler         = [[]];
RealEuler     = [[]];
time          = [];

% init values of omega
Omega_P = [0 0 0];
Omega_I = [0 0 0];

% just to know if we are in the first iteration
first = 1;

tic
i = 1;

fprintf(serialPort, 'o');

% run the loop for 15 seconds
while toc<15
    
    % read a packet
    packet = serialPort.fscanf;
    disp(packet);
    
    % if the packet is not empty
    if ~(isempty(packet))
        
        % try to decode it
        try
            aux = strsplit(packet, '=');
            aux = aux{2};
            fields = strsplit(aux, ',');
            
            % get the time delta
            time(i) = str2double(fields{1});
            
            % get al sensor data
            Accelerometer(i,1) = str2double(fields{2});
            Accelerometer(i,2) = str2double(fields{3});
            Accelerometer(i,3) = str2double(fields{4});
            
            Gyroscope(i,1) = str2double(fields{5});
            Gyroscope(i,2) = str2double(fields{6});
            Gyroscope(i,3) = str2double(fields{7});
            
            Magnetometer(i,1) = str2double(fields{8});
            Magnetometer(i,2) = str2double(fields{9});
            Magnetometer(i,3) = str2double(fields{10});
            
            % get real euler angles calculated by the Arduino
            RealEuler(i,1) = str2double(fields{12}); %pitch
            RealEuler(i,2) = str2double(fields{13}); %roll
            RealEuler(i,3) = str2double(fields{11}); %yaw
            
            % in the first iteration init the dcm matrix
            if (first == 1)
                first = 0;
                
                % init the dcm matrix
                [ pitch, roll, yaw, dcm_init ] = ...
                    reset_fusion( Accelerometer(i,1), ...
                                  Accelerometer(i,2), ...
                                  Accelerometer(i,3), ...
                                  Magnetometer(i,1),  ...
                                  Magnetometer(i,2),  ...
                                  Magnetometer(i,3));
            
            % for the rest of the iterations
            else
                
                % update the dcm matrix with the previous one and the
                % values
                [ pitch, roll, yaw, ...
                  dcm_init, ...
                  Omega_I,  ...
                  Omega_P ] = dcm_algorithm ( dcm_init, ...
                                              Accelerometer(i,1), ...
                                              Accelerometer(i,2), ...
                                              Accelerometer(i,3), ...
                                              Gyroscope(i,1)*(pi/180), ...
                                              Gyroscope(i,2)*(pi/180), ...
                                              Gyroscope(i,3)*(pi/180), ...
                                              Magnetometer(i,1), ...
                                              Magnetometer(i,2), ...
                                              Magnetometer(i,3), ...
                                              time(i), ...
                                              pitch, roll, yaw, ...
                                              Omega_I, Omega_P );
            end
            
            % store the values
            Euler(i,1) = pitch *(180/pi);
            Euler(i,2) = roll  *(180/pi);
            Euler(i,3) = yaw   *(180/pi);
            
            i = i + 1;
        catch error
            disp('Lost packet/Error')
            disp(getReport(error))
        end
    end
    
    % dont care if the packet is empty
    if isempty(packet)
        continue;
    end
    
end

% close the serial port
fclose(serialPort);
delete(serialPort);

% calculate the time axis
time = cumsum(time) - time;

% plot the data
figure('Name', 'Sensor Data');
axis(1) = subplot(3,1,1);
hold on;
plot(time, Gyroscope(:,1), 'r');
plot(time, Gyroscope(:,2), 'g');
plot(time, Gyroscope(:,3), 'b');
legend('X', 'Y', 'Z');
xlabel('Time (s)');
ylabel('Angular rate (deg/s)');
title('Gyroscope');
hold off;
axis(2) = subplot(3,1,2);
hold on;
plot(time, Accelerometer(:,1), 'r');
plot(time, Accelerometer(:,2), 'g');
plot(time, Accelerometer(:,3), 'b');
legend('X', 'Y', 'Z');
xlabel('Time (s)');
ylabel('Acceleration (g)');
title('Accelerometer');
hold off;
axis(3) = subplot(3,1,3);
hold on;
plot(time, Magnetometer(:,1), 'r');
plot(time, Magnetometer(:,2), 'g');
plot(time, Magnetometer(:,3), 'b');
legend('X', 'Y', 'Z');
xlabel('Time (s)');
ylabel('Flux (G)');
title('Magnetometer');
hold off;
linkaxes(axis, 'x');

figure('Name', 'Euler');
hold on
plot(time, Euler(:,1),'r')
plot(time, Euler(:,2),'g')
plot(time, Euler(:,3),'b')
title('Euler Angles Calculated Matlab');
xlabel('Time (s)');
ylabel('Angle (deg)');
legend('Pitch','Roll','Yaw')

figure('Name', 'EulerReal');
hold on
plot(time, RealEuler(:,1),'r')
plot(time, RealEuler(:,2),'g')
plot(time, RealEuler(:,3),'b')
title('Euler Angles Real DCM');
xlabel('Time (s)');
ylabel('Angle (deg)');
legend('Pitch','Roll','Yaw')

% pitch
figure
subplot(311)
hold on
plot(time, Euler(:,1))
plot(time, RealEuler(:,1))
grid on; grid minor;
title('Pitch');
xlabel('Time (s)');
ylabel('Angle (deg)');
legend('Pitch Matlab', 'Pitch Arduino')

% roll
subplot(312)
hold on
plot(time, Euler(:,2))
plot(time, RealEuler(:,2))
grid on; grid minor;
title('Roll');
xlabel('Time (s)');
ylabel('Angle (deg)');
legend('Roll Matlab', 'Roll Arduino')

% yaw
subplot(313)
hold on
plot(time, Euler(:,3))
plot(time, RealEuler(:,3))
grid on; grid minor;
title('Yaw');
xlabel('Time (s)');
ylabel('Angle (deg)');
legend('Yaw Matlab', 'Yaw Arduino')