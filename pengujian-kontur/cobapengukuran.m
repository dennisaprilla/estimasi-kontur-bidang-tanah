close all; clear; clc;

addpath('/home/dennis/Documents/MATLAB/pengujian-kontur/pengukuran1');
addpath('/home/dennis/Documents/MATLAB/pengujian-kontur/pengukuran2');
addpath('/home/dennis/Documents/MATLAB/pengujian-kontur/pengukuran3');

density_gt = 0.1;
pengujian = 1;
jumlah_algoritma = 3;

%% Persiapan GT dan Estimasi

% % Kodingan original, jangan diapus
% x_gt = 0:density_gt:10;
% y_gt = zeros(1, length(x_gt));
% coordinate_gt = [x_gt; y_gt]';
% 
% density_est = 1;
% x_est = 0:density_est:10;
% x_est = x_est + normrnd( 0, 0.2, [1 length(x_est)] )
% y_est = ones(1,length(x_est))
% y_est = y_est + normrnd( 0, 0.1, [1 length(x_est)] );
% coordinate_est = [x_est; y_est]';

if (pengujian==1)
    
    load('pengukuran1_dcm.mat');
    load('pengukuran1_madgwick.mat');
    load('pengukuran1_compFilt.mat');
    load('groundtruth_1.mat');
    
    coordinate_gt       = coords;
    surface_gt          = coords_surfacexyz;
    rollcoordinate_gt   = coords_koordinatroll;
    rollstringvalue_gt  = coords_nilairoll;
    
    coordinate_est1      = pengukuran1_dcm_xyz;
    surface_est1         = pengukuran1_dcm_surfacexyz;
    rollcoordinate_est1  = pengukuran1_dcm_koordinatroll;
    rollstringvalue_est1 = pengukuran1_dcm_nilairoll;
    
    coordinate_est2      = pengukuran1_madgwick_xyz;
    surface_est2         = pengukuran1_madgwick_surfacexyz;
    rollcoordinate_est2  = pengukuran1_madgwick_koordinatroll;
    rollstringvalue_est2 = pengukuran1_madgwick_nilairoll;
    
    coordinate_est3      = pengukuran1_compFilt_xyz;
    surface_est3         = pengukuran1_compFilt_surfacexyz;
    rollcoordinate_est3  = pengukuran1_compFilt_koordinatroll;
    rollstringvalue_est3 = pengukuran1_compFilt_nilairoll;

elseif (pengujian==2)
    
    load('pengukuran2_dcm.mat');
    load('pengukuran2_madgwick.mat');
    load('pengukuran2_compFilt.mat');
    load('groundtruth_2.mat');

    coordinate_gt       = coords;
    surface_gt          = coords_surfacexyz;
    rollcoordinate_gt   = coords_koordinatroll;
    rollstringvalue_gt  = coords_nilairoll;
    
    coordinate_est1      = pengukuran2_dcm_xyz;
    surface_est1         = pengukuran2_dcm_surfacexyz;
    rollcoordinate_est1  = pengukuran2_dcm_koordinatroll;
    rollstringvalue_est1 = pengukuran2_dcm_nilairoll;
    
    coordinate_est2      = pengukuran2_madgwick_xyz;
    surface_est2         = pengukuran2_madgwick_surfacexyz;
    rollcoordinate_est2  = pengukuran2_madgwick_koordinatroll;
    rollstringvalue_est2 = pengukuran2_madgwick_nilairoll;
    
    coordinate_est3      = pengukuran2_compFilt_xyz;
    surface_est3         = pengukuran2_compFilt_surfacexyz;
    rollcoordinate_est3  = pengukuran2_compFilt_koordinatroll;
    rollstringvalue_est3 = pengukuran2_compFilt_nilairoll;
else
    
    load('groundtruth_3.mat');
    load('pengukuran3_dcm.mat');
    load('pengukuran3_madgwick.mat');
    load('pengukuran3_compFilt.mat');

    coordinate_gt       = coords;
    surface_gt          = coords_surfacexyz;
    rollcoordinate_gt   = coords_koordinatroll;
    rollstringvalue_gt  = coords_nilairoll;
    
    coordinate_est1      = pengukuran3_dcm_xyz;
    surface_est1         = pengukuran3_dcm_surfacexyz;
    rollcoordinate_est1  = pengukuran3_dcm_koordinatroll;
    rollstringvalue_est1 = pengukuran3_dcm_nilairoll;
    
    coordinate_est2      = pengukuran3_madgwick_xyz;
    surface_est2         = pengukuran3_madgwick_surfacexyz;
    rollcoordinate_est2  = pengukuran3_madgwick_koordinatroll;
    rollstringvalue_est2 = pengukuran3_madgwick_nilairoll;
    
    coordinate_est3      = pengukuran3_compFilt_xyz;
    surface_est3         = pengukuran3_compFilt_surfacexyz;
    rollcoordinate_est3  = pengukuran3_compFilt_koordinatroll;
    rollstringvalue_est3 = pengukuran3_compFilt_nilairoll;
    
end


%% Cari nearest neighbor GT ke Estimasi buat Perbandingan Kontur Tanah
idx1 = knnsearch(coordinate_gt, coordinate_est1);
idx2 = knnsearch(coordinate_gt, coordinate_est2);
idx3 = knnsearch(coordinate_gt, coordinate_est3);

coordinate_gt2est_1 = coordinate_gt(idx1,:);
coordinate_gt2est_2 = coordinate_gt(idx2,:);
coordinate_gt2est_3 = coordinate_gt(idx3,:);

% % Kodingan original, jangan diapus
% figure(1);
% plot(x_gt, y_gt, '.r');
% grid on; hold on; axis equal;
% plot(x_est, y_est, '.b');
% plot(coordinate_gt2est(:,1), coordinate_gt2est(:,2), 'og');

% % Gambar yang menunjukkan sampel terdekat titik dengan pengukurannya
% figure(1);
% plot3(coordinate_gt(:,1), coordinate_gt(:,2), coordinate_gt(:,3), '.b');
% grid on; hold on; axis equal;
% plot3(coordinate_est1(:,1), coordinate_est1(:,2), coordinate_est1(:,3), '.r');
% plot3(coordinate_gt2est(:,1), coordinate_gt2est(:,2), coordinate_gt2est(:,3), 'og', 'LineWidth', 2);
% title('Sampel Pengukuran dan Titik Terdekat Groundtruthnya');
% legend('Groundtruth', 'Estimation', 'Titik Terdekat');
% xlabel('X (cm)');
% ylabel('Y (cm)');% zlabel('Z (cm)');

figure('Name', 'Evaluasi Kuantitatif: Kontur Tanah', 'unit', 'normalized', 'outerposition', [0 0 1 1]);
title('Perbandingan Estimasi Kontur Tanah Antar Algoritma');
plot3(coordinate_gt(:,1), coordinate_gt(:,2), coordinate_gt(:,3), '-b', 'LineWidth', 3);
grid on; hold on; axis equal;
plot3(coordinate_est1(:,1), coordinate_est1(:,2), coordinate_est1(:,3), '-r', 'LineWidth', 3);
plot3(coordinate_est2(:,1), coordinate_est2(:,2), coordinate_est2(:,3), '-m', 'LineWidth', 3);
plot3(coordinate_est3(:,1), coordinate_est3(:,2), coordinate_est3(:,3), '-g', 'LineWidth', 3);
title('Perbandingan Algoritma Estimasi Kontur Tanah');
legend({'Groundtruth', 'Estimasi DCM-3DC', 'Estimasi Madgwick', '(Baseline) Comp. Filter'}, 'Location', 'best');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');

% view(90,0);
% view(0,90);

%% Visualisasi Perbandingan Bidang Tanah

figure('Name', 'Evaluasi Kuantitatif: Bidang Tanah','units', 'normalized','outerposition',[0 0 1 1]);

surf( surface_gt(:,:,1), surface_gt(:,:,2), surface_gt(:,:,3), ...
      'FaceAlpha',0.5, 'EdgeColor', 'interp'); axis equal; grid on; hold on;
h1 = plot3(coordinate_gt(:,1), coordinate_gt(:,2), coordinate_gt(:,3), '--k'); 

surf( surface_est1(:,:,1), surface_est1(:,:,2), surface_est1(:,:,3), ...
      'FaceAlpha',0.5, 'EdgeColor', 'r');
h2 = plot3(surface_est1(:,:,1), surface_est1(:,:,2), surface_est1(:,:,3), '.k');
h3 = plot3(coordinate_est1(:,1), coordinate_est1(:,2), coordinate_est1(:,3), '-r', 'LineWidth', 3);

surf( surface_est2(:,:,1), surface_est2(:,:,2), surface_est2(:,:,3), ...
      'FaceAlpha',0.5, 'EdgeColor', 'm');
h4 = plot3(surface_est2(:,:,1), surface_est2(:,:,2), surface_est2(:,:,3), '.k');
h5 = plot3(coordinate_est2(:,1), coordinate_est2(:,2), coordinate_est2(:,3), '-m', 'LineWidth', 3);

surf( surface_est3(:,:,1), surface_est3(:,:,2), surface_est3(:,:,3), ...
      'FaceAlpha',0.5, 'EdgeColor', 'g');
h6 = plot3(surface_est3(:,:,1), surface_est3(:,:,2), surface_est3(:,:,3), '.k');
h7 = plot3(coordinate_est3(:,1), coordinate_est3(:,2), coordinate_est3(:,3), '-g', 'LineWidth', 3);
% text(rollcoordinate_est2(:,1), rollcoordinate_est2(:,2), rollcoordinate_est2(:,3), rollstringvalue_est2, ...
%     'VerticalAlignment', 'bottom', ...
%     'HorizontalAlignment', 'center', ...
%     'Position', [0 0 10]);

legend([h1, h3, h5, h7], {'Groundtruth', 'Estimasi DCM-3DC', 'Estimasi Madgwick', '(Baseline) Comp. Filter'}, 'Location', 'best');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');

c = colorbar;
c.Label.String = 'Elevasi relatif terhadap titik awal (cm)';
title('Perbandingan Estimasi Bidang Tanah Antar Algoritma');

% view(90,0);
% view(0,90);

%% Evaluasi Per-Algoritma

% DCM vs Groundtruth ---------------------

figure('Name', 'Evaluasi Kuantitatif: Bidang Tanah DCM-3DC','units', 'normalized','outerposition',[0 0 1 1]);

surf( surface_gt(:,:,1), surface_gt(:,:,2), surface_gt(:,:,3), ...
      'FaceAlpha',0.5, 'EdgeColor', 'interp'); axis equal; grid on; hold on;
h1 = plot3(coordinate_gt(:,1), coordinate_gt(:,2), coordinate_gt(:,3), '--k');

surf( surface_est1(:,:,1), surface_est1(:,:,2), surface_est1(:,:,3), ...
      'FaceAlpha',0.5, 'EdgeColor', 'r');
h2 = plot3(surface_est1(:,:,1), surface_est1(:,:,2), surface_est1(:,:,3), '.k');
h3 = plot3(coordinate_est1(:,1), coordinate_est1(:,2), coordinate_est1(:,3), '-r', 'LineWidth', 3);
text(rollcoordinate_est1(:,1), rollcoordinate_est1(:,2), rollcoordinate_est1(:,3), rollstringvalue_est1, ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'center', ...
    'Position', [0 0 10]);

legend([h1, h3], {'Groundtruth', 'Estimasi DCM-3DC'}, 'Location', 'best');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');

c = colorbar;
c.Label.String = 'Elevasi relatif terhadap titik awal (cm)';
title('Perbandingan Estimasi Bidang Tanah DCM-3DC vs Groundtruth');

% view(90,0);
% view(0,90);


% Madgwick vs Groundtruth -----------------

figure('Name', 'Evaluasi Kuantitatif: Bidang Tanah Madgwick','units', 'normalized','outerposition',[0 0 1 1]);

surf( surface_gt(:,:,1), surface_gt(:,:,2), surface_gt(:,:,3), ...
      'FaceAlpha',0.5, 'EdgeColor', 'interp'); axis equal; grid on; hold on;
h1 = plot3(coordinate_gt(:,1), coordinate_gt(:,2), coordinate_gt(:,3), '--k');

surf( surface_est2(:,:,1), surface_est2(:,:,2), surface_est2(:,:,3), ...
      'FaceAlpha',0.5, 'EdgeColor', 'm');
h4 = plot3(surface_est2(:,:,1), surface_est2(:,:,2), surface_est2(:,:,3), '.k');
h5 = plot3(coordinate_est2(:,1), coordinate_est2(:,2), coordinate_est2(:,3), '-m', 'LineWidth', 3);
text(rollcoordinate_est2(:,1), rollcoordinate_est2(:,2), rollcoordinate_est2(:,3), rollstringvalue_est2, ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'center', ...
    'Position', [0 0 10]);

legend([h1, h5], {'Groundtruth', 'Estimasi Madgwick'}, 'Location', 'best');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');

c = colorbar;
c.Label.String = 'Elevasi relatif terhadap titik awal (cm)';
title('Perbandingan Estimasi Bidang Tanah Madgwick vs Groundtruth');

% view(90,0);
% view(0,90);


% Madgwick vs Groundtruth -----------------

figure('Name', 'Evaluasi Kuantitatif: Bidang Tanah Comp. Filter','units', 'normalized','outerposition',[0 0 1 1]);

surf( surface_gt(:,:,1), surface_gt(:,:,2), surface_gt(:,:,3), ...
      'FaceAlpha',0.5, 'EdgeColor', 'interp'); axis equal; grid on; hold on;
h1 = plot3(coordinate_gt(:,1), coordinate_gt(:,2), coordinate_gt(:,3), '--k');

surf( surface_est3(:,:,1), surface_est3(:,:,2), surface_est3(:,:,3), ...
      'FaceAlpha',0.5, 'EdgeColor', 'g');
h4 = plot3(surface_est3(:,:,1), surface_est3(:,:,2), surface_est3(:,:,3), '.k');
h5 = plot3(coordinate_est3(:,1), coordinate_est3(:,2), coordinate_est3(:,3), '-g', 'LineWidth', 3);
text(rollcoordinate_est3(:,1), rollcoordinate_est3(:,2), rollcoordinate_est3(:,3), rollstringvalue_est3, ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'center', ...
    'Position', [0 0 10]);

legend([h1, h5], {'Groundtruth', '(Baseline) Comp. Filter'}, 'Location', 'best');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');

c = colorbar;
c.Label.String = 'Elevasi relatif terhadap titik awal (cm)';
title('Perbandingan Estimasi Bidang Tanah Comp. Filter vs Groundtruth');

% view(90,0);
% view(0,90);

%% cari total error

RSE(1,:) = sqrt(sum(((coordinate_gt2est_1-coordinate_est1).^2),2));
RSE(2,:) = sqrt(sum(((coordinate_gt2est_2-coordinate_est2).^2),2));
RSE(3,:) = sqrt(sum(((coordinate_gt2est_3-coordinate_est3).^2),2));

for i=1:jumlah_algoritma
    metrik_perbandingan(i,1) = mean(RSE(i,:));
    metrik_perbandingan(i,2) = median(RSE(i,:));
    metrik_perbandingan(i,3) = mean(diff(RSE(i,:)));
    metrik_perbandingan(i,4) = median(diff(RSE(i,:)));
end

% Visualisasi error --------------------------

figure('Name', 'Evaluasi Kuantitatif', 'unit', 'normalized', 'outerposition', [0 0 1 1]);
subplot(2,2,1);
title('Error Growth, Mean, dan Median, DCM-3DC');
hold on; grid on;
plot(RSE(1,:), '-r', 'LineWidth', 2);
plot(ones(1, length(RSE(1,:))) * metrik_perbandingan(1,1), '--b');
plot(ones(1, length(RSE(1,:))) * metrik_perbandingan(1,2), '--g');
legend('Error growth', 'Mean', 'Median');
xlabel('Sampel');
ylabel('Root Square Error (cm)');
ylim([0 max(max(RSE))+5]);

subplot(2,2,2);
title('Error Growth, Mean, dan Median, Madgwick');
hold on; grid on;
plot(RSE(2,:), '-m', 'LineWidth', 2);
plot(ones(1, length(RSE(2,:))) * metrik_perbandingan(2,1), '--b');
plot(ones(1, length(RSE(2,:))) * metrik_perbandingan(2,2), '--g');
legend('Error growth', 'Mean', 'Median');
xlabel('Sampel');
ylabel('Root Square Error (cm)');
ylim([0 max(max(RSE))+5]);

subplot(2,2,3);
title('Error Growth, Mean, dan Median, Comp. Filter');
hold on; grid on;
plot(RSE(3,:), '-g', 'LineWidth', 2);
plot(ones(1, length(RSE(3,:))) * metrik_perbandingan(3,1), '--b');
plot(ones(1, length(RSE(3,:))) * metrik_perbandingan(3,2), '--g');
legend('Error growth', 'Mean', 'Median');
xlabel('Sampel');
ylabel('Root Square Error (cm)');
ylim([0 max(max(RSE))+5]);

%%
% RSE_dcm = sqrt(sum(((coordinate_gt2est_1-coordinate_est1).^2),2));
% RSE_madgwick = sqrt(sum(((coordinate_gt2est_2-coordinate_est2).^2),2));
% RSE_compFilt = sqrt(sum(((coordinate_gt2est_3-coordinate_est3).^2),2));
% 
% distance_RMSE_dcm   = mean(RSE_dcm)
% graph_RMSE_dcm      = ones(1, length(RSE_dcm)) * distance_RMSE_dcm;
% distance_RMedSE_dcm = median(RSE_dcm)
% graph_RMedSE_dcm    = ones(1, length(RSE_dcm)) * distance_RMedSE_dcm;
% 
% distance_RMSE_madgwick   = mean(RSE_madgwick)
% graph_RMSE_madgwick      = ones(1, length(RSE_madgwick)) * distance_RMSE_madgwick;
% distance_RMedSE_madgwick = median(RSE_madgwick)
% graph_RMedSE_madgwick    = ones(1, length(RSE_madgwick)) * distance_RMedSE_madgwick;
% 
% distance_RMSE_compFilt   = mean(RSE_compFilt)
% graph_RMSE_compFilt      = ones(1, length(RSE_compFilt)) * distance_RMSE_compFilt;
% distance_RMedSE_compFilt = median(RSE_compFilt)
% graph_RMedSE_compFilt    = ones(1, length(RSE_compFilt)) * distance_RMedSE_compFilt;

% figure('Name', 'Evaluasi Kuantitatif', 'unit', 'normalized', 'outerposition', [0 0 1 1]);
% subplot(2,2,1);
% title('Error Growth, Mean, dan Median, DCM-3DC');
% hold on; grid on;
% plot(RSE_dcm, '-r', 'LineWidth', 2);
% plot(graph_RMSE_dcm, '--b');
% plot(graph_RMedSE_dcm, '--g');
% legend('Error growth', 'Mean', 'Median');
% xlabel('Sampel');
% ylabel('Root Square Error (cm)');
% ylim([0 max([RSE_dcm; RSE_madgwick; RSE_compFilt])+5]);
% 
% subplot(2,2,2);
% title('Error Growth, Mean, dan Median, Madgwick');
% hold on; grid on;
% plot(RSE_madgwick, '-m', 'LineWidth', 2);
% plot(graph_RMSE_madgwick, '--b');
% plot(graph_RMedSE_madgwick, '--g');
% legend('Error growth', 'Mean', 'Median');
% xlabel('Sampel');
% ylabel('Root Square Error (cm)');
% ylim([0 max([RSE_dcm; RSE_madgwick; RSE_compFilt])+5]);
% 
% subplot(2,2,3);
% title('Error Growth, Mean, dan Median, Comp. Filter');
% hold on; grid on;
% plot(RSE_compFilt, '-g', 'LineWidth', 2);
% plot(graph_RMSE_compFilt, '--b');
% plot(graph_RMedSE_compFilt, '--g');
% legend('Error growth', 'Mean', 'Median');
% xlabel('Sampel');
% ylabel('Root Square Error (cm)');
% ylim([0 max([RSE_dcm; RSE_madgwick; RSE_compFilt])+5]);
