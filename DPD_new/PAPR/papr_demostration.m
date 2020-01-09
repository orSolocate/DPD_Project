close all;
clear;
clc;

%% Initialization
load('100MHzLTE.mat');
load('mask.mat');
x = (waveform(1:50000))';
Smax = sqrt(98/42);
Th = 0.08 * EVM(x, Smax);
Amax_factor = 5
mask_factor = 5
mask = mask_factor*mask;
Amax = Amax_factor*mean(abs(x));

%% Before CFR
figure(1);
plot(abs(x));
plot(abs(fftshift(fft(x))));
hold on;
% plot(mask);
title('FFT of Original Signal and Original Mask');
grid;
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','GridAlpha'),'GridAlpha',1);
papr_x_dB = calcPAPR(x)

%% CFR
x_CFR = PAPR(x, Th, mask, Amax, Smax);

%% After CFR
figure(2);
plot(abs(fftshift(fft(x_CFR))));
hold on;
plot(mask);
title('FFT of Post-CFR Signal and Original Mask');
grid;
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','GridAlpha'),'GridAlpha',1);
papr_xCFR_dB = calcPAPR(x_CFR)

fprintf('Initial PAPR: %fdB, After CFR PAPR: %fdB. Improvment of %fdB\n', papr_x_dB, papr_xCFR_dB, papr_x_dB - papr_xCFR_dB);