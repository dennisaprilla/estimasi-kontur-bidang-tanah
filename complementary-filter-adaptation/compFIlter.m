%% Program Estimasi Kontur Bidang Tanah DCM
% kalo mau simpen nilai ke .mat, jangan lupa ganti nama variablenya
% tergantung pengukuran ke berapa.

clear all; clc; close all;

% pathnya harus disesuaikan
addpath('/home/dennis/Documents/MATLAB/estimasi-kontur-bidang-tanah/complementary-filter-adaptation/quaternion_library');
addpath('/home/dennis/Documents/MATLAB/estimasi-kontur-bidang-tanah/complementary-filter-adaptation/data_pengukuran');

pengukuran = 2;

%%

if (pengukuran==1)
    % teras rumah 3 yang terbaik
    load('groundtruth_1.mat');
    fileID = fopen('teras3.txt','r');

elseif (pengukuran==2)
    % tanjakan bank soal 6 yang terbaik
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

% Magnetometer = A(9:11, :);
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

% Accelerometer = A(3:5,:)' ./ 256;
% Gyroscope = A(6:8,:)' ./ 14.375;
acc = A(2:4,:)' ./ 256;
gyr = A(5:7,:)' ./ 14.375;

time = [1:length(acc)]';

%%

x = 1:length(gyr);

% % groundtruth_2
% gyr_smooth(:,1) = smooth(x, Gyroscope(:,1), 0.1,'sgolay');
% gyr_smooth(:,2) = smooth(x, Gyroscope(:,2), 0.1,'sgolay');
% gyr_smooth(:,3) = smooth(x, Gyroscope(:,3), 0.1,'sgolay');


gyr_smooth(:,1) = smooth(x, gyr(:,1), 0.07,'loess');
gyr_smooth(:,2) = smooth(x, gyr(:,2), 0.07,'loess');
gyr_smooth(:,3) = smooth(x, gyr(:,3), 0.07,'loess');

% % groundtruth_2s
% acc_smooth(:,1) = smooth(x, Accelerometer(:,1), 0.025,'sgolay');
% acc_smooth(:,2) = smooth(x, Accelerometer(:,2), 0.025,'sgolay');
% acc_smooth(:,3) = smooth(x, Accelerometer(:,3), 0.025,'sgolay');

acc_smooth(:,1) = smooth(x, acc(:,1), 0.03,'moving');
acc_smooth(:,2) = smooth(x, acc(:,2), 0.03,'moving');
acc_smooth(:,3) = smooth(x, acc(:,3), 0.03,'moving');

% % groundtruth_2
% mag_smooth(:,1) = smooth(x, Magnetometer(:,1), 0.05,'sgolay');
% mag_smooth(:,2) = smooth(x, Magnetometer(:,2), 0.05,'sgolay');
% mag_smooth(:,3) = smooth(x, Magnetometer(:,3), 0.05,'sgolay');

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


% subplot(3,1,2);
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
ylabel('Flux (Gs)');
title('Magnetometer');
hold off;

% %%
% figure('Name', 'Sensor Data');
% subplot(3,1,1);
% hold on;
% plot(Gyroscope(:,1), 'r');
% plot(Gyroscope(:,2), 'g');
% plot(Gyroscope(:,3), 'b');
% legend('X', 'Y', 'Z');
% xlabel('Sampel');
% ylabel('Kecepatan Sudut (deg/s)');
% title('Gyroscope');
% hold off;
% 
% subplot(3,1,2);
% hold on;
% plot(Accelerometer(:,1), 'r');
% plot(Accelerometer(:,2), 'g');
% plot(Accelerometer(:,3), 'b');
% legend('X', 'Y', 'Z');
% xlabel('Sampel');
% ylabel('Akselerasi (g)');
% title('Accelerometer');
% hold off;
% 
% subplot(3,1,3);
% hold on;
% plot(Magnetometer(:,1), 'r');
% plot(Magnetometer(:,2), 'g');
% plot(Magnetometer(:,3), 'b');
% legend('X', 'Y', 'Z');
% xlabel('Sampel');
% ylabel('Flux (Gs)');
% title('Magnetometer');
% hold off;

%%

% replace variable
gyr=gyr_smooth;
acc=acc_smooth;
mag=mag_smooth;

% inisialisasi
roll = 0;
pitch = 0;
yaw = 0;

% jeda waktu rata2
dt = 1/121;
tau = 0.98;

eulers = [];

for i=1:length(gyr)

    % % Find angles from accelerometer
    % accelPitch = rad2deg(atan2(acc(i,2), acc(i,3)));
    % accelRoll = rad2deg(atan2(acc(i,1), acc(i,3)));
    
    accelRoll = rad2deg(atan2(acc(i,2), acc(i,3)));
    accelPitch = rad2deg(atan2(acc(i,1), acc(i,3)));
    
    % Apply complementary filter
    pitch = (tau)*(pitch + gyr(i,1) * dt) + (1 - tau)*(accelPitch);
    roll = (tau)*(roll - gyr(i,2) * dt) + (1 - tau)*(accelRoll);
    yaw = (yaw + gyr(i,3) * dt);
    
    % eulers = [eulers; roll, pitch, yaw];
    eulers = [eulers; roll, pitch, -yaw];
end

eulers = deg2rad(eulers);

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

%% script tambahan

% % groundtruth_1 & groundtruth_2
% d_jarak = A(2,:)' .* 1.97;
d_jarak = A(1,:)' .* 1.94;

% konversi matrix jadi euler
% eulers = rotMat2euler(R, 'yxz');

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
% karena dcm make orientasi global (bumi), yaw harus diubah jadi lokal (sensor)

% proyeksi ke sumbu X dan Y
projected_dx = sin(yaws).*projected_djarak;
projected_dy = cos(yaws).*projected_djarak;

% trajektori ground, terdiri dari X dan Y
trajektori_x = cumsum(projected_dx);
trajektori_y = cumsum(projected_dy);

xyz = [trajektori_x, trajektori_y , trajektori_altitude];

% display
figure('Name', 'Trajectory', 'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.r');
axis equal; grid on; hold on;
plot3(coords(:,1), coords(:,2), coords(:,3), '.b');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
legend('Estimasi','Groundtruth')
title('Estimasi Trayektori / Kontur Tanah');
view(0,0);

subplot(1,2,2);
plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.r');
axis equal; grid on; hold on;
plot3(coords(:,1), coords(:,2), coords(:,3), '.b');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
legend('Estimasi','Groundtruth')
title('Estimasi Trayektori / Kontur Tanah');
view(90,90);

%% script estimasi bidang tanah

resolusi_bidangtanah = 200;
% resolusi_bidangtanah = 100;

rolls = eulers(:,1);
roll_eulers = [-yaws rolls zeros(length(rolls),1)];
roll_matrices = eul2rotm(roll_eulers);

point_for_surface = repmat([-10 0 0; 10 0 0]', 1, 1, length(gyr));
translation = reshape(xyz', 3, 1, []);

surfaces = [];

for i=1:resolusi_bidangtanah:length(gyr)    
    a = roll_matrices(:,:,i)*point_for_surface(:,:,i)+translation(:,:,i);
    surfaces = [surfaces; a'];
end

% surfaces_unique = unique(surfaces, 'rows');
x_surface = reshape(surfaces(:,1), 2, []);
y_surface = reshape(surfaces(:,2), 2, []);
z_surface = reshape(surfaces(:,3), 2, []);

% tampilin
figure('Name', 'Hasil Trajektori');
s = surf(x_surface,y_surface,z_surface,'FaceAlpha',0.8);
axis equal; grid on; hold on;

plot3(x_surface, y_surface, z_surface, '.k');
plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.r');

% coba coba untuk bisa ngasih teks nilai roll ke dalem plot 3D
data = 1:resolusi_bidangtanah:length(R);
koordinat_data = xyz(data,:);
nilairoll_data = round(rad2deg(eulers(data,1)), 1);
nilairoll_datastring = cellstr(num2str(nilairoll_data));
text(koordinat_data(:,1), koordinat_data(:,2), koordinat_data(:,3), nilairoll_datastring, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Position', [0 0 10]);

% rename variable buat disimpen ke .mat (buat perbandingan)
pengukuran3_compFilt_xyz = xyz;
pengukuran3_compFilt_surfacexyz(:,:,1) = x_surface;
pengukuran3_compFilt_surfacexyz(:,:,2) = y_surface;
pengukuran3_compFilt_surfacexyz(:,:,3) = z_surface;
pengukuran3_compFilt_koordinatroll = koordinat_data;
pengukuran3_compFilt_nilairoll = nilairoll_datastring;

axis equal;

c = colorbar;
c.Label.String = 'Elevasi relatif terhadap titik awal (cm)';

xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
title('Trayektori + Estimasi Surface');

% save('pengukuran3_compFilt.mat', ...
%      'pengukuran3_compFilt_xyz', ...
%      'pengukuran3_compFilt_surfacexyz', ...
%      'pengukuran3_compFilt_koordinatroll', ...
%      'pengukuran3_compFilt_nilairoll');