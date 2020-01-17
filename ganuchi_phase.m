clc; clear all; close all
%% Comment are for matlab amplifier
%amp_data          = 'simrfV2_powamp_dpd_data.mat';
%load(amp_data);                                            %load amp data
signal            = 'y_WL_100MHzLTE.mat';
load(signal);                                              %load signal
x=waveform(start_pos_sig:start_pos_sig+X_length-1); %if you load the file:
%y_WL_100MHzLTE.mat
%x                     = x(1:5000)./(2*norm(x,2));
%RMS_in                = 10*log10( norm(x)^2/50/length(x)) + 30 + 20;
%[y, ~, ~, ~] = RFWebLab_PA_meas_v1_1(x, RMS_in);
%mem_depth = 4 ;                                        
%mem_deg   = 5 ;  
%AMP_coef_Matrix  = Get_coef_MP(inDataPA, outDataPA, mem_deg, mem_depth);
%y     = [zeros(mem_depth,1); Get_model_output_MP(AMP_coef_Matrix, x, mem_deg, mem_depth)];
%%
% WL_delay = finddelay(x, y);
% if(WL_delay >= 0)
%     y = y(WL_delay+1:end);
%     x = x(1:end-WL_delay);
% elseif(WL_delay < 0)
%     y = y(1:end+WL_delay);
%     x = x(1-WL_delay:end);
% end
%%
N = length(x);
% window preparation
win = rectwin(N);
% fft of the first signal
X = fft(x.*win);
% fft of the second signal
Y = fft(y.*win);

%fs=SR

Fv = linspace(0, 1, fix(length(x)/2)+1)*fs/2;           % Frequency Vector
Iv = 1:length(Fv);

phs_x = angle(X);
phsv_x = phs_x(Iv); 
phs_y = angle(Y);
phsv_y = phs_y(Iv); 

phi=log(phsv_y./phsv_x);
freq=Fv/fs;
freq=freq';

dphi=diff(phi)./diff(freq);

figure();
plot (freq,phi);
%hold on;
%plot(freq(2:end),dphi);
title('phase/freq')

P = polyfit(freq,phi,1);
slope = P(1)
yFit = polyval(P, freq);
hold on;
plot(freq, yFit, 'b', 'LineWidth', 2);
legend('phi(f/fs)','linear fit');


t0=-1*real(slope)/(2*pi)

%T=1/fs 
%t0=t0_normalized/fs

Samples_delay=round(t0)