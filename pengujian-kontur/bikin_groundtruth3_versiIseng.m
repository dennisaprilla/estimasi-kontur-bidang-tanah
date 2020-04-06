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

%%%

step2a_x = 0 : -density : (-86)+(-650)+(-83);
step2a_y = ones(1, length(step2a_x)) * (220);
step2a_z = zeros(1, length(step2a_x));

step3a_y = (220) : density : (220)+(166);
step3a_x = ones(1, length(step3a_y)) * ((-86)+(-650)+(-83));
step3a_z = zeros(1, length(step3a_y));

step4a_x = (-86)+(-650)+(-83) : density : 0;
step4a_y = ones(1, length(step4a_x)) * ((220)+(166));
step4a_z = zeros(1, length(step4a_x));

%%%

step4b_x = (-86)+(-650)+(-83) : density : 0;
step4b_y = ones(1, length(step4b_x)) * ((220)+(166));
step4b_z = ones(1, length(step4b_x)) * (80);

%%%

coord1 = [step1_x; step1_y; step1_z];
coord2 = [step2_x; step2_y; step2_z];
coord3 = [step3_x; step3_y; step3_z];
coord4 = [step4_x; step4_y; step4_z];
coord5 = [step5_x; step5_y; step5_z];
coord6 = [step6_x; step6_y; step6_z];
coord7 = [step7_x; step7_y; step7_z];
coord8 = [step8_x; step8_y; step8_z];

coords2a = [step2a_x; step2a_y; step2a_z];
coords3a = [step3a_x; step3a_y; step3a_z];
coords4a = [step4a_x; step4a_y; step4a_z];


%%%

coords = [coord1, coord2, coord3, coord4, coord5, coord6, coord7, coord8]';

coordsa = [coord1, coords2a, coords3a, coords4a]';
coordsb = [step4b_x; step4b_y; step4b_z]';

figure(1);
plot3(coordsa(:,1), coordsa(:,2), coordsa(:,3), '-k'); axis equal; hold on;
plot3(coordsb(:,1), coordsb(:,2), coordsb(:,3), '-k');
plot3(coords(:,1), coords(:,2), coords(:,3), '.b');
title('Groundtruth Kontur Tanah');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
