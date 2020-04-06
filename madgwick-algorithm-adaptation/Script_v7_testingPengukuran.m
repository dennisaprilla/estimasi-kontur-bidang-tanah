%% Housekeeping
 
addpath('ximu_matlab_library');	% include x-IMU MATLAB library
addpath('quaternion_library');	% include quatenrion library
addpath('data_pengukuran');	    % include quatenrion library
close all;                     	% close all figures
clear;                         	% clear all variables
clc;                          	% clear the command terminal

pengukuran = 1;

%% Import data

if (pengukuran==1)
    % teras rumah 3 yang terbaik
    load('groundtruth_1.mat');
    fileID = fopen('teras3.txt','r');

elseif (pengukuran==2)
    % tanjakan bank soal 5 yang terbaik
    load('groundtruth_2.mat');
    fileID = fopen('tanjakan_banksoal_pelan6.txt','r');

else
    % tanjakan luar 5 yang terbaik
    load('groundtruth_3.mat');
    fileID = fopen('005.TXT','r');
end

% formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f\n';
% sizeA = [11 Inf];
% A = fscanf(fileID,formatSpec,sizeA);
% fclose(fileID);

formatSpec = '%f %f %f %f %f %f %f %f %f %f %f\n';
sizeA = [10 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);

% mag = A(9:11, :);
mag = A(8:10, :);

if (pengukuran==1)
    % rumah pak imam
    t = [-1420.204027; -632.604186; -443.972665];
    R = [0.798837, -0.063148, 0.019997; ...
        -0.063148, 0.784889, -0.018052; ...
        0.019997, -0.018052, 0.922340];
else 
    % bagian soal
    t = [-1274.824572; -695.609169; -146.653054];
    R = [0.920123, -0.068939, 0.016445; ...
        -0.068939, 0.895958, -0.021982; ...
        0.016445, -0.021982, 1.061931];
end

Magnetometer_translated = bsxfun(@minus, mag, t);
Magnetometer_calibrated = R*Magnetometer_translated;
mag = mag' ./ 12000;

% acc = A(3:5,:)' ./ 256;
% gyr = A(6:8,:)' ./ 14.375;
acc = A(2:4,:)' ./ 256;
gyr = A(5:7,:)' ./ 14.375;

time = [1:length(acc)]';

%%
x = 1:length(gyr);

% % groundtruth_2
% gyr_smooth(:,1) = smooth(x, gyr(:,1), 0.1,'sgolay');
% gyr_smooth(:,2) = smooth(x, gyr(:,2), 0.1,'sgolay');
% gyr_smooth(:,3) = smooth(x, gyr(:,3), 0.1,'sgolay');

gyr_smooth(:,1) = smooth(x, gyr(:,1), 0.07,'loess');
gyr_smooth(:,2) = smooth(x, gyr(:,2), 0.07,'loess');
gyr_smooth(:,3) = smooth(x, gyr(:,3), 0.07,'loess');

% % groundtruth_2s
% acc_smooth(:,1) = smooth(x, acc(:,1), 0.025,'sgolay');
% acc_smooth(:,2) = smooth(x, acc(:,2), 0.025,'sgolay');
% acc_smooth(:,3) = smooth(x, acc(:,3), 0.025,'sgolay');

acc_smooth(:,1) = smooth(x, acc(:,1), 0.03,'moving');
acc_smooth(:,2) = smooth(x, acc(:,2), 0.03,'moving');
acc_smooth(:,3) = smooth(x, acc(:,3), 0.03,'moving');

% % groundtruth_2
% mag_smooth(:,1) = smooth(x, mag(:,1), 0.05,'sgolay');
% mag_smooth(:,2) = smooth(x, mag(:,2), 0.05,'sgolay');
% mag_smooth(:,3) = smooth(x, mag(:,3), 0.05,'sgolay');

mag_smooth(:,1) = smooth(x, mag(:,1), 0.03,'moving');
mag_smooth(:,2) = smooth(x, mag(:,2), 0.03,'moving');
mag_smooth(:,3) = smooth(x, mag(:,3), 0.03,'moving');


%figure('Name', 'Sensor Data');
% subplot(3,1,1);
figure(1);
hold on;
plot(gyr(:,1), 'r');
plot(gyr_smooth(:,1), 'c', 'LineWidth', 2);
plot(gyr(:,2), 'g');
plot(gyr_smooth(:,2), 'm', 'LineWidth', 2);
plot(gyr(:,3), 'b');
plot(gyr_smooth(:,3), 'y', 'LineWidth', 2);
grid on;
legend('X', 'Xs', 'Y', 'Ys', 'Z', 'Zs');
xlabel('Sampel');
ylabel('Kecepatan Sudut (deg/s)');
title('Gyroscope');
hold off;

% subplot(3,1,2);
figure(2);
hold on;
plot(acc(:,1), 'r');
plot(acc_smooth(:,1), 'c', 'LineWidth', 2);
plot(acc(:,2), 'g');
plot(acc_smooth(:,2), 'm', 'LineWidth', 2);
plot(acc(:,3), 'b');
plot(acc_smooth(:,3), 'y', 'LineWidth', 2);
grid on;
legend('X', 'Xs', 'Y', 'Ys', 'Z', 'Zs');
xlabel('Sampel');
ylabel('Akselerasi (g)');
title('Accelerometer');
hold off;

figure(3);
hold on;
plot(mag(:,1), 'r');
plot(mag_smooth(:,1), 'c', 'LineWidth', 2);
plot(mag(:,2), 'g');
plot(mag_smooth(:,2), 'm', 'LineWidth', 2);
plot(mag(:,3), 'b');
plot(mag_smooth(:,3), 'y', 'LineWidth', 2);
grid on;
legend('X', 'Xs', 'Y', 'Ys', 'Z', 'Zs');
xlabel('Sampel');
ylabel('Akselerasi (g)');
title('Accelerometer');
hold off;

% subplot(3,1,3);
% hold on;
% plot(mag(:,1), 'r');
% plot(mag_smooth(:,1), 'c', 'LineWidth', 2);
% plot(mag(:,2), 'g');
% plot(mag_smooth(:,2), 'm', 'LineWidth', 2);
% plot(mag(:,3), 'b');
% plot(mag_smooth(:,3), 'y', 'LineWidth', 2);
% legend('X', 'Xs', 'Y', 'Ys', 'Z', 'Zs');
% grid on;
% xlabel('Sampel');
% ylabel('Flux (Gs)');
% title('Magnetometer');
% hold off;

% figure('Name', 'Sensor Data');
% subplot(3,1,1);
% hold on;
% plot(mag(:,1), 'r');
% plot(mag_smooth1(:,1), 'b');
% grid on;
% legend('Normal', 'loess', 'sgolay');
% xlabel('Sampel');
% ylabel('Kecepatan Sudut (deg/s)');
% title('Gyroscope X');
% hold off;
% 
% subplot(3,1,2);
% hold on;
% plot(mag(:,2), 'r');
% plot(mag_smooth2(:,2), 'b');
% grid on;
% legend('Normal', 'loess', 'sgolay');
% xlabel('Sampel');
% ylabel('Kecepatan Sudut (deg/s)');
% title('Gyroscope Y');
% hold off;
% 
% subplot(3,1,3);
% hold on;
% plot(mag(:,3), 'r');
% plot(mag_smooth3(:,3), 'b');
% legend('Normal', 'loess', 'sgolay');
% grid on;
% xlabel('Sampel');
% ylabel('Kecepatan Sudut (deg/s)');
% title('Gyroscope Z');
% hold off;


%% Process data through AHRS algorithm (calcualte orientation)
% See: http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/

gyr=gyr_smooth;
acc=acc_smooth;
mag=mag_smooth;

% samplePeriode tergantung dari seberapa cepet arduino nyimpen ngeproses
% dan nyimpen data ke SD card

% ini yang seharusnya, rata2 8-9ms
% samplePeriod = 1/123;
samplePeriod = 1/121;

% ini yang terbagus
% AHRS = MadgwickAHRS('SamplePeriod', samplePeriod, 'Beta', 0.01);

AHRS = MadgwickAHRS('SamplePeriod', samplePeriod, 'Beta', 0.0001);

% buat nampung rotation matrixnya
R = zeros(3,3,length(gyr));

for i = 1:length(gyr)
    AHRS.Update(gyr(i,:) * (pi/180), acc(i,:), mag(i,:));
    
    %pake ini yang terbagus
    %AHRS.UpdateIMU(gyr(i,:) * (pi/180), acc(i,:));
    
    % transpose because ahrs provides Earth relative to sensor
    R(:,:,i) = quatern2rotMat(AHRS.Quaternion)';
end

eulers = rad2deg(rotMat2euler(R));

figure('Name', 'Euler Angles');
hold on; grid on;
plot(time, eulers(:,1), 'r');
plot(time, eulers(:,2), 'g');
plot(time, eulers(:,3), 'b');
title('Euler angles');
xlabel('Time (s)');
ylabel('Angle (deg)');
legend('\phi', '\theta', '\psi');
hold off;

%% script estimasi kontur tanah

% bikin data-dataan nilai jarak antar sampling (satuannya belum ada)
% konstanta tergantung rotary dan ban dari alat
% d_jarak = A(2,:)' .* 1.97;
d_jarak = A(1,:)' .* 1.94;

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

xyz = [-trajektori_x, trajektori_y , -trajektori_altitude];

% rename variable buat disimpen ke .mat (buat perbandingan)
pengukuran1_madgwick_xyz = xyz;

% figure(3);
% subplot(1,2,1);
% plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.r');
% axis equal; grid on; hold on;
% view(90,0);
% 
% subplot(1,2,2);
% plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.r');
% axis equal; grid on; hold on;
% view(90,90);
% 
% figure(4)
% plot(-pitched_deg, '-r');
% grid on;

figure('Name', 'Trajectory','units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.r');
axis equal; grid on; hold on;
plot3(coords(:,1), coords(:,2), coords(:,3), '.b');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% legend('Estimasi','Groundtruth')
title('Estimasi Trayektori / Kontur Tanah');
view(0,0);

subplot(1,2,2);
plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.r');
axis equal; grid on; hold on;
plot3(coords(:,1), coords(:,2), coords(:,3), '.b');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
% legend('Estimasi','Groundtruth')
title('Estimasi Trayektori / Kontur Tanah');
view(90,90);

%% script estimasi bidang tanah

rolls = eulers(:,1);
% roll_eulers = [-yaws rolls zeros(length(rolls),1)];
roll_eulers = [yaws rolls zeros(length(rolls),1)];
roll_matrices = eul2rotm(roll_eulers);

point_for_surface = repmat([-10 0 0; 10 0 0]', 1, 1, length(R));
translation = reshape(xyz', 3, 1, []);

surfaces = [];

for i=1:200:length(R)    
    a = roll_matrices(:,:,i)*point_for_surface(:,:,i)+translation(:,:,i);
    surfaces = [surfaces; a'];
end

% surfaces_unique = unique(surfaces, 'rows');
x_surface = reshape(surfaces(:,1), 2, []);
y_surface = reshape(surfaces(:,2), 2, []);
z_surface = reshape(surfaces(:,3), 2, []);

% rename variable buat disimpen ke .mat (buat perbandingan)
pengukuran1_madgwick_surfacexyz(:,:,1) = x_surface;
pengukuran1_madgwick_surfacexyz(:,:,2) = y_surface;
pengukuran1_madgwick_surfacexyz(:,:,3) = z_surface;

figure('Name', 'Hasil Trajektori', 'units', 'normalized','outerposition',[0 0 1 1]);
s = surf(x_surface,y_surface,z_surface,'FaceAlpha',0.8);
axis equal; grid on; hold on;

plot3(x_surface, y_surface, z_surface, '.k');
plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.r');

% coba coba untuk bisa ngasih teks nilai roll ke dalem plot 3D
data = 1:200:length(R);
koordinat_data = xyz(data,:);
nilairoll_data = round(rad2deg(eulers(data,1)), 1);
nilairoll_datastring = cellstr(num2str(nilairoll_data));
text(koordinat_data(:,1), koordinat_data(:,2), koordinat_data(:,3), nilairoll_datastring, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Position', [0 0 10]);

% rename variable buat disimpen ke .mat (buat perbandingan)
pengukuran1_madgwick_koordinatroll = koordinat_data;
pengukuran1_madgwick_nilairoll = nilairoll_datastring;

axis equal;

c = colorbar;
c.Label.String = 'Elevasi relatif terhadap titik awal (cm)';

xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
title('Trayektori + Estimasi Surface');

% save('pengukuran1_madgwick.mat', 'pengukuran1_madgwick_xyz', 'pengukuran1_madgwick_surfacexyz', 'pengukuran1_madgwick_koordinatroll', 'pengukuran1_madgwick_nilairoll');

%% Play animation

% SamplePlotFreq = 8;
% 
% linPosHP = zeros(length(gyr),3);
% 
% SixDOFanimation(linPosHP, R, ...
%                 'SamplePlotFreq', SamplePlotFreq, 'Trail', 'DotsOnly', ...
%                 'FullScreen', true, ...
%                 'Position', [9 39 1280 720], ...
%                 'AxisLength', 2, 'ShowArrowHead', true, ...
%                 'Xlabel', 'X (m)', 'Ylabel', 'Y (m)', 'Zlabel', 'Z (m)', 'ShowLegend', false, 'Title', 'Unfiltered',...
%                 'CreateAVI', false, 'AVIfileNameEnum', false, 'AVIfps', ((1/samplePeriod) / SamplePlotFreq));            