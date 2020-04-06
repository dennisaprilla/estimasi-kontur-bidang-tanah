%% Program Estimasi Kontur Bidang Tanah DCM
% kalo mau simpen nilai ke .mat, jangan lupa ganti nama variablenya
% tergantung pengukuran ke berapa.

clear all; clc; close all;

% ini pathnya harus diganti ke direktori ini
addpath('/home/dennis/Documents/MATLAB/estimasi-kontur-bidang-tanah/dcm-algorithm-adaptation/src');
addpath('/home/dennis/Documents/MATLAB/estimasi-kontur-bidang-tanah/dcm-algorithm-adaptation/quaternion_library');
addpath('/home/dennis/Documents/MATLAB/estimasi-kontur-bidang-tanah/dcm-algorithm-adaptation/data_pengukuran');

pengukuran = 2;

%%

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

% Magnetometer = A(9:11, :);
Magnetometer = A(8:10, :);

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

Magnetometer_translated = bsxfun(@minus, Magnetometer, t);
Magnetometer_calibrated = R*Magnetometer_translated;
Magnetometer = Magnetometer' ./ 12000;

% Accelerometer = A(3:5,:)' ./ 256;
% Gyroscope = A(6:8,:)' ./ 14.375;
Accelerometer = A(2:4,:)' ./ 256;
Gyroscope = A(5:7,:)' ./ 14.375;

time = [1:length(Accelerometer)]';

%%
x = 1:length(Gyroscope);

% % groundtruth_2
% gyr_smooth(:,1) = smooth(x, Gyroscope(:,1), 0.1,'sgolay');
% gyr_smooth(:,2) = smooth(x, Gyroscope(:,2), 0.1,'sgolay');
% gyr_smooth(:,3) = smooth(x, Gyroscope(:,3), 0.1,'sgolay');


gyr_smooth(:,1) = smooth(x, Gyroscope(:,1), 0.07,'loess');
gyr_smooth(:,2) = smooth(x, Gyroscope(:,2), 0.07,'loess');
gyr_smooth(:,3) = smooth(x, Gyroscope(:,3), 0.07,'loess');

% % groundtruth_2s
% acc_smooth(:,1) = smooth(x, Accelerometer(:,1), 0.025,'sgolay');
% acc_smooth(:,2) = smooth(x, Accelerometer(:,2), 0.025,'sgolay');
% acc_smooth(:,3) = smooth(x, Accelerometer(:,3), 0.025,'sgolay');

acc_smooth(:,1) = smooth(x, Accelerometer(:,1), 0.03,'moving');
acc_smooth(:,2) = smooth(x, Accelerometer(:,2), 0.03,'moving');
acc_smooth(:,3) = smooth(x, Accelerometer(:,3), 0.03,'moving');

% % groundtruth_2
% mag_smooth(:,1) = smooth(x, Magnetometer(:,1), 0.05,'sgolay');
% mag_smooth(:,2) = smooth(x, Magnetometer(:,2), 0.05,'sgolay');
% mag_smooth(:,3) = smooth(x, Magnetometer(:,3), 0.05,'sgolay');

mag_smooth(:,1) = smooth(x, Magnetometer(:,1), 0.03,'moving');
mag_smooth(:,2) = smooth(x, Magnetometer(:,2), 0.03,'moving');
mag_smooth(:,3) = smooth(x, Magnetometer(:,3), 0.03,'moving');

%figure('Name', 'Sensor Data');
% subplot(3,1,1);
figure(1);
hold on;
plot(Gyroscope(:,1), 'r');
plot(gyr_smooth(:,1), 'c', 'LineWidth', 2);
plot(Gyroscope(:,2), 'g');
plot(gyr_smooth(:,2), 'm', 'LineWidth', 2);
plot(Gyroscope(:,3), 'b');
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
plot(Accelerometer(:,1), 'r');
plot(acc_smooth(:,1), 'c', 'LineWidth', 2);
plot(Accelerometer(:,2), 'g');
plot(acc_smooth(:,2), 'm', 'LineWidth', 2);
plot(Accelerometer(:,3), 'b');
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
plot(Magnetometer(:,1), 'r');
plot(mag_smooth(:,1), 'c', 'LineWidth', 2);
plot(Magnetometer(:,2), 'g');
plot(mag_smooth(:,2), 'm', 'LineWidth', 2);
plot(Magnetometer(:,3), 'b');
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
Gyroscope=gyr_smooth;
Accelerometer=acc_smooth;
Magnetometer=mag_smooth;

% buat nampung R matrix
R = zeros(3,3,length(Gyroscope));

% init values of omega
Omega_P = [0 0 0];
Omega_I = [0 0 0];

yaw_init = 0;

for i = 1:length(Gyroscope)
    if i==1
        % init the dcm matrix
        [ pitch, roll, yaw, dcm_init ] = reset_fusion( Accelerometer(i,1), ...
                                                       Accelerometer(i,2), ...
                                                       Accelerometer(i,3), ...
                                                       Magnetometer(i,1),  ...
                                                       Magnetometer(i,2),  ...
                                                       Magnetometer(i,3));
        
        yaw_init = yaw;
        
        rot_roll = eul2rotm([0, 0, roll]);
        rot_pitch = eul2rotm([0, pitch, 0]);
        rot_rollpitch = rot_roll*rot_pitch;
        dcm_init = rot_rollpitch;
        
    else
        
        % update the dcm matrix with the previous one and the values
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
                                      1/121, ...
                                      pitch, roll, yaw, ...
                                      Omega_I, Omega_P );
    end
    
    R(:,:,i)=dcm_init';
end

eulers = rad2deg(rotMat2euler(R, 'yxz'));

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

%% Play animation

% SamplePlotFreq = 8;
% 
% linPosHP = zeros(length(Gyroscope),3);
% 
% SixDOFanimation(linPosHP, R, ...
%                 'SamplePlotFreq', SamplePlotFreq, 'Trail', 'DotsOnly', ...
%                 'Position', [9 39 1280 720], ...
%                 'AxisLength', 0.1, 'ShowArrowHead', false, ...
%                 'Xlabel', 'X (m)', 'Ylabel', 'Y (m)', 'Zlabel', 'Z (m)', 'ShowLegend', false, 'Title', 'Unfiltered',...
%                 'CreateAVI', false, 'AVIfileNameEnum', false, 'AVIfps', ((1/166) / SamplePlotFreq));            

%% script tambahan

% % groundtruth_1 & groundtruth_2
% d_jarak = A(2,:)' .* 1.97;
d_jarak = A(1,:)' .* 1.94;

% konversi matrix jadi euler
eulers = rotMat2euler(R, 'yxz');

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

% rename variable buat disimpen ke .mat (buat perbandingan)
pengukuran2_dcm_xyz = xyz;

% figure('Name', 'Trajectory');
% subplot(1,2,1);
% plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.r');
% axis equal; grid on; hold on;
% view(90,0);
% 
% subplot(1,2,2);
% plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.r');
% axis equal; grid on; hold on;
% view(90,90);


figure('Name', 'Trajectory', 'units','normalized','outerposition',[0 0 1 1]);
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

resolusi_bidangtanah = 200;

rolls = eulers(:,1);
roll_eulers = [-yaws rolls zeros(length(rolls),1)];
roll_matrices = eul2rotm(roll_eulers);

point_for_surface = repmat([-10 0 0; 10 0 0]', 1, 1, length(R));
translation = reshape(xyz', 3, 1, []);

surfaces = [];

for i=1:resolusi_bidangtanah:length(R)    
    a = roll_matrices(:,:,i)*point_for_surface(:,:,i)+translation(:,:,i);
    surfaces = [surfaces; a'];
end

% surfaces_unique = unique(surfaces, 'rows');
x_surface = reshape(surfaces(:,1), 2, []);
y_surface = reshape(surfaces(:,2), 2, []);
z_surface = reshape(surfaces(:,3), 2, []);

% rename variable buat disimpen ke .mat (buat perbandingan)
pengukuran2_dcm_surfacexyz(:,:,1) = x_surface;
pengukuran2_dcm_surfacexyz(:,:,2) = y_surface;
pengukuran2_dcm_surfacexyz(:,:,3) = z_surface;

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
pengukuran2_dcm_koordinatroll = koordinat_data;
pengukuran2_dcm_nilairoll = nilairoll_datastring;

axis equal;

c = colorbar;
c.Label.String = 'Elevasi relatif terhadap titik awal (cm)';

xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
title('Trayektori + Estimasi Surface');

% save('pengukuran2_dcm.mat', 'pengukuran2_dcm_xyz', 'pengukuran2_dcm_surfacexyz', 'pengukuran2_dcm_koordinatroll', 'pengukuran2_dcm_nilairoll');