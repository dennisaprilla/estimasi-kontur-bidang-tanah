close all; clear; clc;

load('DCMErrors.mat');
load('madgwickErrors.mat');
load('DCMDuration.mat');
load('madgwickDuration.mat');

%%

figure('Name', 'Errors');
grid on; hold on;
plot(noise_levels, madgwickErrors_meanObservation, '-or');
plot(noise_levels, DCMErrors_meanObservation, '-^b');
legend('Madgwick', 'DCM');
xlabel('Noise terhadap STD axis');
ylabel('Error dalam Quaternions');
title('Perbandingan Error');

%%

string = {'Madwick'; 'DCM'};

figure('Name', 'Durations');
grid on; hold on;

b = bar([madgwickDurationTimes_mean; DCMDurationTimes_mean]);
set(gca, 'XTickLabel', string, 'XTick',1:numel(string));

l = {'Murni', '+ Penyimpanan', '+ Transformasi'};
legend(b, l);

xlabel('Algoritma');
ylabel('Waktu (detik)');
title('Perbandingan Durasi');

%%

madgwickfile = whos('-file','madgwickDuration.mat');
DCMfile = whos('-file','DCMDuration.mat');

figure('Name', 'Memory Size');
grid on; hold on;
bar([madgwickfile(3).bytes; DCMfile(3).bytes]);
set(gca, 'XTickLabel', string, 'XTick',1:numel(string))
xlabel('Algoritma');
ylabel('Besar Memory yang Dibutuhkan (bytes)');
title('Perbandingan Besar Memory');