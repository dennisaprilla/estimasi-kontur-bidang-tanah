%% Housekeeping

addpath('ximu_matlab_library');	% include x-IMU MATLAB library
addpath('quaternion_library');	% include quatenrion library
close all;                     	% close all figures
clear;                         	% clear all variables
clc;                          	% clear the command terminal

%% Import data

% xIMUdata = xIMUdataClass('LoggedData/LoggedData');
% xIMUdata = xIMUdataClass('Datasets/stairsAndCorridor');
xIMUdata = xIMUdataClass('Datasets/spiralStairs');

samplePeriod = 1/256;

gyr = [xIMUdata.CalInertialAndMagneticData.Gyroscope.X...
       xIMUdata.CalInertialAndMagneticData.Gyroscope.Y...
       xIMUdata.CalInertialAndMagneticData.Gyroscope.Z];        % gyroscope
acc = [xIMUdata.CalInertialAndMagneticData.Accelerometer.X...
       xIMUdata.CalInertialAndMagneticData.Accelerometer.Y...
       xIMUdata.CalInertialAndMagneticData.Accelerometer.Z];	% accelerometer
mag = [xIMUdata.CalInertialAndMagneticData.Magnetometer.X...
       xIMUdata.CalInertialAndMagneticData.Magnetometer.Y...
       xIMUdata.CalInertialAndMagneticData.Magnetometer.Z];     % magnetometer

clear('xIMUdata');


mag = mag - mean(mag);

figure('Name', 'Sensor Data');
subplot(3,1,1);
hold on;
plot(gyr(:,1), 'r');
plot(gyr(:,2), 'g');
plot(gyr(:,3), 'b');
legend('X', 'Y', 'Z');
xlabel('Sampel');
ylabel('Kecepatan Sudut (deg/s)');
title('Gyroscope');
hold off;

subplot(3,1,2);
hold on;
plot(acc(:,1), 'r');
plot(acc(:,2), 'g');
plot(acc(:,3), 'b');
legend('X', 'Y', 'Z');
xlabel('Sampel');
ylabel('Akselerasi (g)');
title('Accelerometer');
hold off;

subplot(3,1,3);
hold on;
plot(mag(:,1), 'r');
plot(mag(:,2), 'g');
plot(mag(:,3), 'b');
legend('X', 'Y', 'Z');
xlabel('Sampel');
ylabel('Flux (Gs)');
title('Magnetometer');
hold off;

   
%% Process data through AHRS algorithm (calcualte orientation)
% See: http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/

R = zeros(3,3,length(gyr));     % rotation matrix describing sensor relative to Earth

AHRS = MadgwickAHRS('SamplePeriod', samplePeriod, 'Beta', 0.0001);
%AHRS = MahonyAHRS('SamplePeriod', samplePeriod, 'Kp', 1);

for i = 1:length(gyr)
    AHRS.Update(gyr(i,:) * (pi/180), acc(i,:), mag(i,:));
    %AHRS.UpdateIMU(gyr(i,:) * (pi/180), acc(i,:));
    R(:,:,i) = quatern2rotMat(AHRS.Quaternion)';    % transpose because ahrs provides Earth relative to sensor
end

%% script tambahan

% bikin data-dataan nilai jarak antar sampling (satuannya belum ada)
d_jarak = sqrt(randn(size(acc,1),1).^2)*0.0025; 

% konversi matrix jadi euler
eulers = rotMat2euler(R);

% ambil pitchnya
pitches = eulers(:,2);
pitched_deg = rad2deg(pitches);

% proyeksi tinggi
projected_dtinggi = sin(pitches).*d_jarak;

% proyeksi ke ground
projected_djarak = cos(pitches).*d_jarak;

% trajektori tinggi terhadap tanah dari pitch
trajektori_altitude = cumsum(projected_dtinggi);

% trajektori di ground yang hanya dari pitch
trajektori_azimuth = cumsum(projected_djarak);

% ambil yaw nya
yaws = eulers(:,3);

% proyeksi ke sumbu X dan Y
projected_dx = sin(yaws).*projected_djarak;
projected_dy = cos(yaws).*projected_djarak;

% trajektori ground, terdiri dari X dan Y
trajektori_x = cumsum(projected_dx);
trajektori_y = cumsum(projected_dy);

xyz = [trajektori_x, trajektori_y , trajektori_altitude];

%% Plot

% figure('Name', 'YZ-Plane');
% hold on; grid on;
% plot(xyz(:,2), xyz(:,3), '-g');
% xlabel('Sumbu Y (meter)');
% ylabel('Sumbu Z (meter)');
% title('Proyeksi Trajektori ke YZ-Plane');
% axis equal;
% 
% [maxvalue, argmax] = max(xyz(:,3));
% yMeter_whenZmax = xyz(argmax,2);
% [minvalue, argmin] = min(xyz(:,3));
% yMeter_whenZmin = xyz(argmin,2);
% 
% plot(yMeter_whenZmax, maxvalue, 'om');
% plot(yMeter_whenZmin, minvalue, 'om');

%%%%%%%%%%%%%%%

% figure('Name', 'XY-Plane');
% hold on; grid on;
% plot(xyz(:,1), xyz(:,2), '-r');
% xlabel('Sumbu X (meter)');
% ylabel('Sumbu Y (meter)');
% title('Proyeksi Trajektori ke XY-Plane');
% axis equal;
% 
% [maxvalue, argmax] = max(xyz(:,1));
% yMeter_whenXmax = xyz(argmax,2);
% [minvalue, argmin] = min(xyz(:,1));
% yMeter_whenXmin = xyz(argmin,2);
% 
% plot(maxvalue, yMeter_whenXmax, 'om');
% plot(minvalue, yMeter_whenXmin, 'om');

%%%%%%%%%%%%%%%

% figure('Name', 'XZ-Plane');
% hold on; grid on;
% plot(xyz(:,1), xyz(:,3), '-b');
% xlabel('Sumbu X (meter)');
% ylabel('Sumbu Z (meter)');
% title('Proyeksi Trajektori ke XZ-Plane');
% axis equal;

% [maxvalue, argmax] = max(xyz(:,1));
% zMeter_whenXmax = xyz(argmax,3);
% [minvalue, argmin] = min(xyz(:,1));
% zMeter_whenXmin = xyz(argmin,3);
% 
% plot(maxvalue, zMeter_whenXmax, 'om');
% plot(minvalue, zMeter_whenXmin, 'om');
% 
% [maxvalue, argmax] = max(xyz(:,3));
% xMeter_whenXmax = xyz(argmax,1);
% [minvalue, argmin] = min(xyz(:,3));
% xMeter_whenXmin = xyz(argmin,1);
% 
% plot(xMeter_whenXmax, maxvalue, '*m');
% plot(xMeter_whenXmin, minvalue, '*m');

%% Plot lagi
jarak_untuk_plot = cumsum(d_jarak);

% figure('Name', 'Perubahan X');
% hold on; grid on;
% plot(jarak_untuk_plot, projected_dx, '-r');
% xlabel('Jarak (meter)');
% ylabel('dX');
% title('Perubahan X');
% 
% figure('Name', 'Perubahan Y');
% hold on; grid on;
% plot(jarak_untuk_plot, projected_dy, '-g');
% xlabel('Jarak (meter)');
% ylabel('dY');
% title('Perubahan Y');
% 
% figure('Name', 'Perubahan Z');
% hold on; grid on;
% plot(jarak_untuk_plot, projected_dtinggi,'-b');
% xlabel('Jarak (meter)');
% ylabel('dZ');
% title('Perubahan Z');

%% Tambahan yang Pak Musa pengen

rolls = eulers(:,1);
%roll_eulers = [-yaws zeros(length(rolls),1) rolls];
roll_eulers = [-yaws rolls zeros(length(rolls),1)];
roll_matrices = eul2rotm(roll_eulers);

point_for_surface = repmat([-0.25 0 0; 0.25 0 0]', 1, 1, length(R));
translation = reshape(xyz', 3, 1, []);

surfaces = [];

for i=1:200:length(R)
%for i=5000:200:9500
%for i=2200:200:10500
%for i=700:200:10300
%for i=1:length(R)
    
    %disp(i);
    %disp(point_for_surface(:,:,i));
    %disp(roll_matrices(:,:,i)*point_for_surface(:,:,i));
    
    a = roll_matrices(:,:,i)*point_for_surface(:,:,i)+translation(:,:,i);
    %surfaces = cat(3,surfaces,a);
    surfaces = [surfaces; a'];
end

% surfaces_unique = unique(surfaces, 'rows');

x_surface = reshape(surfaces(:,1), 2, []);
y_surface = reshape(surfaces(:,2), 2, []);
z_surface = reshape(surfaces(:,3), 2, []);

figure('Name', 'Hasil Trajektori');
s = surf(x_surface,y_surface,z_surface,'FaceAlpha',0.8);
axis equal; grid on; hold on;

plot3(x_surface, y_surface, z_surface, '.k');
plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.r');
%plot3(xyz(5000:9500, 1), xyz(5000:9500, 2), xyz(5000:9500, 3), '.r');
%plot3(xyz(2200:10500, 1), xyz(2200:10500, 2), xyz(2200:10500, 3), '.r');
%plot3(xyz(700:10300, 1), xyz(700:10300, 2), xyz(700:10300, 3), '.r');

% coba coba untuk bisa ngasih teks nilai roll ke dalem plot 3D
data = 1:400:length(R);
%data = 5000:200:9500;
%data = 2200:200:10500;
%data = 700:200:10300;
koordinat_data = xyz(data,:);
nilairoll_data = round(rad2deg(eulers(data,1)), 1);
nilairoll_datastring = cellstr(num2str(nilairoll_data));
text(koordinat_data(:,1), koordinat_data(:,2), koordinat_data(:,3), nilairoll_datastring, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Position', [0 0 10]);

axis equal;
view(52, 30);
%view(90,90);
%view(90,0);
%view(0,0);
%view(120, 40);

c = colorbar;
c.Label.String = 'Elevasi relatif terhadap titik awal (meter)';

xlabel('X (meter)');
ylabel('Y (meter)');
zlabel('Z (meter)');
title('Trayektori + Estimasi Surface');