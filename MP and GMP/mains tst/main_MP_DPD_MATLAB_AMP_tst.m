clc; clear; close all;                                     %start freash

%% General User inputs
signal            = '100MHzLTE.mat';
model_run_period  = 60000 ;                                %num of samples to run in the model
start_pos_sig     = 1;
end_pos_sig       = start_pos_sig+model_run_period-1;

%% Loads
load(signal);                                              %load signal
load('inDataPA.mat');
load('outDataPA.mat');

%% Model User inputs
iterations_num_MP  = 20;
resistor           = 50;
miu_MP             = -0.1;

mem_depth = 2 ;                                            %M in the MP model
mem_deg   = 5 ;                                            %K in the MP model

%% Initialize and Model
avg_gain         = abs(mean(outDataPA./inDataPA));         
x                = waveform(start_pos_sig:end_pos_sig);
y_d              = x.*avg_gain;
AMP_coef_Matrix  = Get_coef_MP(inDataPA, outDataPA, mem_deg, mem_depth);
PD_coef_Matrix   = Get_coef_MP(outDataPA./avg_gain, inDataPA, mem_deg, mem_depth);

%% initial calculations
x_opt_MP  = [zeros(mem_depth,1); PD_MP(x./avg_gain, PD_coef_Matrix, mem_deg, mem_depth)];

y_MP     = [zeros(mem_depth,1); Get_model_output_MP(AMP_coef_Matrix, x_opt_MP, mem_deg, mem_depth)];
y_MP     = y_MP.*avg_gain;
y_in_MP  = [zeros(mem_depth,1); Get_model_output_MP(AMP_coef_Matrix, x, mem_deg, mem_depth)];

%phase correction
%x_opt_MP  = ifft(fft(x_opt_MP).*exp(-phdiffmeasure(y_MP, x_opt_MP)*1i));

%% Sub calculations
transfer_no_PD        = y_in_MP./x;
transfer_with_MP_PD   = y_MP./x;

%% Ploting initial state
figure
pspectrum(y_d(mem_depth+1:end))
hold on
pspectrum(y_in_MP(mem_depth+1:end))
hold on
pspectrum(y_MP(mem_depth+1:end))
legend('y_d','no PD in MP','with MP PD')

figure;
plot(abs(x), 20*log10(abs(transfer_no_PD)), 'o');
hold on;
plot(abs(x_opt_MP), 20*log10(abs(transfer_with_MP_PD)), '*');
legend('Without DPD in MP','with MP PD')
xlabel('abs(input)');
ylabel('Gain');
title('Linearity Comparison');

%% Iterate pre
error_per_iter_MP         = ones(iterations_num_MP+1,1).*inf;

phase_per_iter_MP         = zeros(iterations_num_MP+1,1);
phase_per_iter_MP(1)      = phdiffmeasure(y_MP, x_opt_MP)*180/pi;

avg_power_per_iter_MP     = zeros(iterations_num_MP+1,1);
avg_power_per_iter_MP(1)  = 10*log10( norm(x_opt_MP)^2/resistor/length(x)) + 30;

PAPR_per_iter_MP          = zeros(iterations_num_MP+1,1);
PAPR_per_iter_MP(1)       = 20*log10(max(abs(x_opt_MP))*sqrt(length(x_opt_MP))/norm(x_opt_MP));

RMS_in_per_iter_MP        = zeros(iterations_num_MP+1,1);
RMS_in_per_iter_MP(1)     = avg_power_per_iter_MP(1);

%% Iterate MP
ll = 1;
while (ll<=iterations_num_MP)
    
    y_MP                    = ifft(fft(y_MP).*exp(-phdiffmeasure(y_d, y_MP)*1i));
    output_err_MP           = y_d - y_MP;
    error_per_iter_MP(ll+1) = norm(output_err_MP,2);
     
    %get the updated coefficients for the model
    PD_coef_Matrix = Update_coef_MP(y_MP./avg_gain, x_opt_MP, PD_coef_Matrix, mem_deg, mem_depth, miu_MP);
    AMP_coef_Matrix  = Update_coef_MP(y_MP, x_opt_MP, AMP_coef_Matrix, mem_deg, mem_depth, miu_MP);
    
    %get the updated optimal input
    x_opt_MP       = [zeros(mem_depth,1); PD_MP(x./avg_gain, PD_coef_Matrix, mem_deg, mem_depth)];
    
    %phase correction
    %x_opt_MP       = ifft(fft(x_opt_MP).*exp(-phdiffmeasure(y_MP, x_opt_MP)*1i));
    
    %get the Amp output
    y_MP  = [zeros(mem_depth,1); Get_model_output_MP(AMP_coef_Matrix, x_opt_MP, mem_deg, mem_depth)];
    y_MP  = y_MP.*avg_gain;
    
    %spectrum plot in each 2nd iteration
    if(0 == mod(ll,2))
        figure
        pspectrum(y_d(mem_depth+1:end))
        hold on
        pspectrum(y_in_MP(mem_depth+1:end))
        hold on
        pspectrum(y_MP(mem_depth+1:end))
        legend('y_d','no PD','with MP PD')
    end
       
    disp(['MP iteration num: ',num2str(ll), ' is done.'])
    ll = ll+1;
    
    phase_per_iter_MP(ll)     = phdiffmeasure(y_MP, x_opt_MP)*180/pi;
    avg_power_per_iter_MP(ll) = 10*log10( norm(x_opt_MP)^2/resistor/length(x)) + 30;
    PAPR_per_iter_MP(ll)      = 20*log10(max(abs(x_opt_MP))*sqrt(length(x_opt_MP))/norm(x_opt_MP));
    RMS_in_per_iter_MP(ll)    = RMS_in_per_iter_MP(ll-1) - (PAPR_per_iter_MP(ll) - PAPR_per_iter_MP(ll-1));
  
end

%% Sub calculations
transfer_no_PD        = y_in_MP./x;
transfer_with_MP_PD   = y_MP./x;

%% Ploting final results
figure;
plot(abs(x), 20*log10(abs(transfer_no_PD)), 'o');
hold on;
plot(abs(x_opt_MP), 20*log10(abs(transfer_with_MP_PD)), '*');
legend('Without DPD in MP','with MP PD')
xlabel('abs(input)');
ylabel('Gain');
title('Linearity Comparison');

figure
plot(error_per_iter_MP(2:end))
legend('error per iteration MP')
xlabel('iteration #');
ylabel('error on y (norm-2)');

figure
plot(phase_per_iter_MP)
legend('phase per iteration MP')
xlabel('iteration #');
ylabel('phase diff [deg]');

figure
plot(avg_power_per_iter_MP)
legend('avg power per iteration MP')
xlabel('iteration #');
ylabel('average input power [dBm]');

figure
plot(RMS_in_per_iter_MP)
legend('RMS in per iteration MP')
xlabel('iteration #');
ylabel('RMS in [dBm]');

figure
pspectrum(y_d(mem_depth+1:end))
hold on
pspectrum(y_in_MP(mem_depth+1:end))
hold on
pspectrum(y_MP(mem_depth+1:end))
legend('y_d','no PD in MP','with MP PD')
