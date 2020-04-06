clear;
close all;
clc;

% ganti density kalo mau lebih rapat lagi
density = 0.1;

%% Bikin Groundtruth Koordinat

step1_y = 0:density:220;
step1_x = zeros(1, length(step1_y));
step1_z = zeros(1, length(step1_y));

step2_x = (0) : -density : (-86);
step2_y = ones(1, length(step2_x)) * (220);
step2_z = zeros(1, length(step2_x));

step3_x = (-86) : -density : (-86)+(-650);
step3_y = ones(1, length(step3_x))*220;
step3_z = linspace(0, 80, length(step3_x));

step4_x = (-86)+(-650) : -density : (-86)+(-650)+(-83);
step4_y = ones(1, length(step4_x)) * (220);
step4_z = ones(1, length(step4_x)) * (80);

step5_y = (220) : density : (220)+(166);
step5_x = ones(1, length(step5_y)) * ((-86)+(-650)+(-83));
step5_z = ones(1, length(step5_y)) * (80);

step6_x = (-86)+(-650)+(-83) : density : (-86)+(-650);
step6_y = ones(1, length(step6_x)) * ((220)+(166));
step6_z = ones(1, length(step6_x)) * (80);

step7_x = (-86)+(-650) : density : (-86);
step7_y = ones(1, length(step7_x)) * ((220)+(166));
step7_z = linspace(80, (80)+(80), length(step7_x));

step8_x = (-86) : density : 0;
step8_y = ones(1, length(step8_x)) * ((220)+(166));
step8_z = ones(1, length(step8_x)) * ((80)+(80));

coord1 = [step1_x; step1_y; step1_z];
coord2 = [step2_x; step2_y; step2_z];
coord3 = [step3_x; step3_y; step3_z];
coord4 = [step4_x; step4_y; step4_z];
coord5 = [step5_x; step5_y; step5_z];
coord6 = [step6_x; step6_y; step6_z];
coord7 = [step7_x; step7_y; step7_z];
coord8 = [step8_x; step8_y; step8_z];

coords = [coord1, coord2, coord3, coord4, coord5, coord6, coord7, coord8]';

x_sampel = 1:length(coords);
coords_smooth(:,1) = smooth(x_sampel, coords(:,1), 0.02,'moving');
coords_smooth(:,2) = smooth(x_sampel, coords(:,2), 0.02,'moving');
coords_smooth(:,3) = smooth(x_sampel, coords(:,3), 0.02,'moving');

%% Bikin Groundtruth Roll-Yaw

rolls = zeros(1,length(coords))';

yaw1 = ones(1, length(coord1)) * 0;
yaw2 = ones(1, length([coord2, coord3, coord4])) * 90;
yaw3 = ones(1, length(coord5)) * 0;
yaw4 = ones(1, length([coord6, coord7, coord8])) * -90;

yaws = [yaw1, yaw2, yaw3, yaw4];
yaws_smooth = deg2rad(smooth(x_sampel, yaws, 0.02,'moving'));

roll_eulers = [ yaws_smooth, rolls, zeros(length(rolls),1) ];
roll_matrices = eul2rotm(roll_eulers, 'zyx');

figure(1);
plot3(coords(:,1), coords(:,2), coords(:,3), '.b');
grid on;
axis equal;
title('Groundtruth Kontur Tanah');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');

%%
resolusi_bidangtanah = 300;

point_for_surface = repmat([-50 0 0; 50 0 0]', 1, 1, length(rolls));
translation = reshape(coords_smooth', 3, 1, []);

surfaces = [];

for i=1:resolusi_bidangtanah:length(rolls)    
    a = roll_matrices(:,:,i)*point_for_surface(:,:,i)+translation(:,:,i);
    surfaces = [surfaces; a'];
end

% surfaces_unique = unique(surfaces, 'rows');
x_surface = reshape(surfaces(:,1), 2, []);
y_surface = reshape(surfaces(:,2), 2, []);
z_surface = reshape(surfaces(:,3), 2, []);

% rename variable buat disimpen ke .mat (buat perbandingan)
coords_surfacexyz(:,:,1) = x_surface;
coords_surfacexyz(:,:,2) = y_surface;
coords_surfacexyz(:,:,3) = z_surface;

figure('Name', 'Hasil Trajektori', 'units', 'normalized','outerposition',[0 0 1 1]);
surf(x_surface,y_surface,z_surface,'FaceAlpha',0.8); axis equal; grid on; hold on;
plot3(x_surface, y_surface, z_surface, '.k');
plot3(coords_smooth(:,1), coords_smooth(:,2), coords_smooth(:,3), '.r');

% coba coba untuk bisa ngasih teks nilai roll ke dalem plot 3D
data = 1:resolusi_bidangtanah:length(rolls);
koordinat_data = coords(data,:);
nilairoll_data = round(rad2deg(rolls(data,1)), 1);
nilairoll_datastring = cellstr(num2str(nilairoll_data));
text(koordinat_data(:,1), koordinat_data(:,2), koordinat_data(:,3), nilairoll_datastring, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Position', [0 0 10]);

% rename variable buat disimpen ke .mat (buat perbandingan)
coords_koordinatroll = koordinat_data;
coords_nilairoll = nilairoll_datastring;

c = colorbar;
c.Label.String = 'Elevasi relatif terhadap titik awal (cm)';

xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
title('Trayektori + Estimasi Surface');

%% iseng

% figure(1);
% subplot(1,2,1);
% plot3(coords(:,1), coords(:,2), coords(:,3), '-b', 'LineWidth', 2); grid on; axis equal;
% subplot(1,2,2);
% plot3(coords_smooth(:,1), coords_smooth(:,2), coords_smooth(:,3), '-', 'LineWidth', 2, 'Color', [0.3010 0.7450 0.9330]); grid on; axis equal;
% title('Groundtruth Kontur Tanah');
% xlabel('X (cm)');
% ylabel('Y (cm)');
% zlabel('Z (cm)');

% figure(2);
% plot(coords(:,2), '-r'); grid on; hold on;
% plot(yaws_smooth, '-b');
% line([length(coord1) length(coord1)], [-100 400]);
% title('Groundtruth Kontur Tanah');