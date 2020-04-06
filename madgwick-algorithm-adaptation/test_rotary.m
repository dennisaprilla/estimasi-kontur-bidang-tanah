%% asumsi:
% kecepatan alatnya rata-rata 50 cm/s
% rotary memberikan output: 360 pulse per revolution
% diameter roda alat: 6 cm

% berapa kecepatan rata-rata kalo satuannya pulse per second?
% 50/37.7 * 360 = 477.5 pulse/s

% kalau rata2 kecepatan 50 cm/s berapa pulse yang terjadi dalam waktu 1 detik
% 1 detik = 50 cm = 477.5 pulse

% kalau rata2 kecepatan 50 cm/s, 1 pulse berapa jaraknya?
% 50 cm / 477.5 pulse = 0.104 cm/pulse = 377/3600

% kalau rata2 kecepatan 50 cm/s, 1 pulse berapa detik?
% 1 s / 477.5 pulse = 0.002094 s/pulse = 377/180000

%% bikin data dummy

clear; close all; clc;

pulse_in_cm = 377/3600;
pulse_in_s  = 377/180000;
n_sample = 10000;

d_jarak = repmat(pulse_in_cm, 1, n_sample);
d_jarak_noise = repmat(pulse_in_cm, 1,n_sample) + normrnd(0, 0.01, [1, n_sample]);

d_t  = repmat(pulse_in_s, 1, n_sample);
time = cumsum(d_t);

%% gambar2

sample_to_plot = 100;
 
figure(1);

% Create the first axes
hax1 = axes();

% Plot something here
hplot1 = stem(d_jarak_noise(1:sample_to_plot), 'Color', 'r');
axis([0 sample_to_plot 0 0.2]);

grid on;

% get handle to current axes
hax1 = gca;
% set box property to off 
set(hax1,'box','off','color','white')

hax2 = axes('Position', get(hax1, 'Position'),'box','off', ...  % Copy position
            'XAxisLocation', 'top', ...             % Put the x axis on top
            'YAxisLocation', 'right', ...           % Doesn't really matter           
            'Color', 'none', ...                    % Make it transparent
            'YTick', []); 


hplot2 = line(time(1:sample_to_plot), d_jarak(1:sample_to_plot), 'Color', 'b', 'Parent', hax2);
axis([0 time(sample_to_plot) 0 0.2]);

% Link the y limits and position together
linkprop([hax1, hax2], {'ylim', 'Position'});

% Draw some labels
xlabel(hax1, 'Blue Line');
xlabel(hax2, 'Red Line');
ylabel(hax1, 'Some Value');

% Add a legend? Why not?!
legend([hplot1, hplot2], {'Perubahan Jarak + Noise', 'Perubahan Jarak Rata-rata'})