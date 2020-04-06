clear;
close all;
clc;

% ganti density kalo mau lebih rapat lagi
density = 0.1;

step1_y = 0:density:200;
step1_x = zeros(1, length(step1_y));
step1_z = zeros(1, length(step1_y));

step2_x = 0:-density:-45;
step2_y = ones(1, length(step2_x))*200;
step2_z = zeros(1, length(step2_x));

step3_x = (-45):-density:(-45)+(-318);
step3_y = ones(1, length(step3_x))*200;
step3_z = linspace(0, 50, length(step3_x));

step4_x = (-45)+(-318):-density:(-45)+(-318)+(-52);
step4_y = ones(1, length(step4_x))*200;
step4_z = ones(1, length(step4_x))*50;

step5_y = 200:-density:0;
step5_x = ones(1, length(step5_y))*( (-45)+(-318)+(-52) );
step5_z = ones(1, length(step5_y))*50;

step2a_x = 0:-density:(-45)+(-318)+(-52);
step2a_y = ones(1, length(step2a_x))*200;
step2a_z = zeros(1, length(step2a_x));

step3a_y = 200:-density:0;
step3a_x = ones(1, length(step3a_y))*( (-45)+(-318)+(-52) );
step3a_z = zeros(1, length(step3a_y));


coord1 = [step1_x; step1_y; step1_z];
coord2 = [step2_x; step2_y; step2_z];
coord3 = [step3_x; step3_y; step3_z];
coord4 = [step4_x; step4_y; step4_z];
coord5 = [step5_x; step5_y; step5_z];

coord2a = [step2a_x; step2a_y; step2a_z];
coord3a = [step3a_x; step3a_y; step3a_z];

coords = [coord1, coord2, coord3, coord4, coord5]';
coordsa = [coord1, coord2a, coord3a]';

figure(1);
plot3(coordsa(:,1), coordsa(:,2), coordsa(:,3), '-k'); hold on; axis equal;
plot3(coords(:,1), coords(:,2), coords(:,3), '.b');
% grid on;
title('Groundtruth Kontur Tanah');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');