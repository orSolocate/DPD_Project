clc; clear; close all;                                     %start freash

%% General User inputs
%signal            = '100MHzLTE.mat';
signal            = 'testsignal.mat';
model_run_period  = 60000 ;                                %num of samples to run in the model
start_pos_sig     = 1;
end_pos_sig       = start_pos_sig+model_run_period-1;
Fs                = 200e6;

%% Model User inputs
iterations_num_MP  = 4;
resistor           = 50;
miu_MP             = 0.5;

mem_depth = 4 ;                                            %M in the MP model
mem_deg   = 5 ;                                            %K in the MP model

%% Loads
load(signal);                                              %load signal
%x                     = waveform(start_pos_sig:end_pos_sig);
x                     = x(start_pos_sig:end_pos_sig)./norm(x,2);
RMS_in                = 10*log10( norm(x)^2/resistor/length(x)) + 30;
[y, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_1(x);

%% Initialize and Model
WL_delay            = finddelay(x, y);
if(WL_delay >=0)
    avg_gain        = abs(mean(y(WL_delay+1:end)./x(1:end-WL_delay)));
    PD_coef_MP_Mat  = Get_coef_MP((y(WL_delay+1:end)')./avg_gain, x(1:end-WL_delay)', mem_deg, mem_depth);
    AMP_coef_Matrix = Get_coef_MP(x(1:end-WL_delay)', y(WL_delay+1:end)', mem_deg, mem_depth);
end
if(WL_delay < 0)
    avg_gain        = abs(mean(y(1:end-WL_delay)./x(WL_delay+1:end)));
    PD_coef_MP_Mat  = Get_coef_MP((y(1:end-WL_delay)')./avg_gain, x(WL_delay+1:end)', mem_deg, mem_depth);
    AMP_coef_Matrix = Get_coef_MP(x(WL_delay+1:end)', y(1:end-WL_delay)', mem_deg, mem_depth);
end
y_d                 = x.*avg_gain;

%% initial calculations
x_opt_MP  = [zeros(mem_depth,1); PD_MP(x./avg_gain, PD_coef_MP_Mat, mem_deg, mem_depth)];
x_opt_MP  = [zeros(mem_depth,1); PD_MP(x_opt_MP, AMP_coef_Matrix, mem_deg, mem_depth)];

[y_MP, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_1(x_opt_MP);
    
%phase correction
%x_opt_MP  = ifft(fft(x_opt_MP).*exp(-phdiffmeasure(y_MP, x_opt_MP)*1i));

%% Sub calculations
WL_delay1 = finddelay(x,y);
transfer_no_PD        = y(WL_delay1+1:end)./x(1:end-WL_delay1);
WL_delay2 = finddelay(x,y_MP);
transfer_with_MP_PD   = y_MP(WL_delay2+1:end)./x(1:end-WL_delay2);

%% Ploting initial state
figure
pspectrum(y_d(mem_depth+1:end), Fs)
hold on
pspectrum(y(mem_depth+1:end), Fs)
hold on
pspectrum(y_MP(mem_depth+1:end), Fs)
legend('y_d','no PD','with MP PD')

figure;
plot(abs(x(1:end-WL_delay1)), 20*log10(abs(transfer_no_PD)), 'o');
hold on;
plot(abs(x_opt_MP(1:end-WL_delay2)), 20*log10(abs(transfer_with_MP_PD)), '*');
legend('Without DPD','with MP PD')
xlabel('abs(input)');
ylabel('Gain');
title('Linearity Comparison');

%% Iterate pre
error_per_iter_MP         = ones(iterations_num_MP+1,1).*inf;

phase_per_iter_MP         = zeros(iterations_num_MP+1,1);
phase_per_iter_MP(1)      = phdiffmeasure(y_MP, x_opt_MP)*180/pi;

avg_power_per_iter_MP     = zeros(iterations_num_MP+1,1);
avg_power_per_iter_MP(1)  = 10*log10( norm(x_opt_MP)^2/resistor/length(x_opt_MP)) + 30;

PAPR_per_iter_MP          = zeros(iterations_num_MP+1,1);
PAPR_per_iter_MP(1)       = 20*log10(max(abs(x_opt_MP))*sqrt(length(x_opt_MP))/norm(x_opt_MP));

RMS_in_per_iter_MP        = zeros(iterations_num_MP+1,1);
RMS_in_per_iter_MP(1)     = avg_power_per_iter_MP(1);

err                       = zeros(iterations_num_MP+1,1);
err(1)                    = mean(abs(y_MP(1:end-WL_delay) - x_opt_MP(WL_delay+1:end)));
%% Iterate MP
ll = 1;
while (ll<=iterations_num_MP)
    
    y_MP                     = ifft(fft(y_MP).*exp(-phdiffmeasure(y_d, y_MP)*1i));
    output_err_MP            = Get_output_err_vec(y_d, y_MP);
    error_per_iter_MP(ll+1)  = norm(output_err_MP,2);
    
    %get the updated coefficients for the model
    PD_coef_MP_Mat  = Update_coef_MP(y_MP./avg_gain, x_opt_MP, PD_coef_MP_Mat, mem_deg, mem_depth, miu_MP);
    AMP_coef_Matrix = Update_coef_MP(x_opt_MP, y_MP, AMP_coef_Matrix, mem_deg, mem_depth, miu_MP);
   
    %get the updated optimal input
    x_opt_MP = [zeros(mem_depth,1); PD_MP(x./avg_gain, PD_coef_MP_Mat, mem_deg, mem_depth)];
    x_opt_MP = [zeros(mem_depth,1); PD_MP(x_opt_MP, AMP_coef_Matrix, mem_deg, mem_depth)];
    
    %phase correction
    %x_opt_MP                 = ifft(fft(x_opt_MP).*exp(-phdiffmeasure(y_MP, x_opt_MP)*1i));
    
    %get the Amp output
    [y_MP, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_1(x_opt_MP);
    
    %spectrum plot every 2nd iteration
    if(0 == mod(ll,1))
        figure
        pspectrum(y_d(mem_depth+1:end), Fs)
        hold on
        pspectrum(y(mem_depth+1:end), Fs)
        hold on
        pspectrum(y_MP(mem_depth+1:end), Fs)
        %hold on
        %pspectrum(x_opt_MP(mem_depth+1:end))
        legend('y_d','no PD','with MP PD', 'x_opt')
    end
       
    disp(['MP iteration num: ',num2str(ll), ' is done.'])
    ll = ll+1;
    
    phase_per_iter_MP(ll)     = phdiffmeasure(y_MP, x_opt_MP)*180/pi;
    avg_power_per_iter_MP(ll) = 10*log10( norm(x_opt_MP)^2/resistor/length(x_opt_MP)) + 30;
    err(ll)                   = mean(abs(y_MP(1:end-WL_delay) - x_opt_MP(WL_delay+1:end)));
    PAPR_per_iter_MP(ll)      = 20*log10(max(abs(x_opt_MP))*sqrt(length(x_opt_MP))/norm(x_opt_MP));
    RMS_in_per_iter_MP(ll)    = RMS_in_per_iter_MP(ll-1) - (PAPR_per_iter_MP(ll) - PAPR_per_iter_MP(ll-1));

end

%% Sub calculations
WL_delay = finddelay(x,y);
transfer_no_PD        = y(WL_delay+1:end)./x(1:end-WL_delay);
WL_delay = finddelay(x,y_MP);
transfer_with_MP_PD   = y_MP(WL_delay+1:end)./x(1:end-WL_delay);

%% Ploting final results
figure;
plot(abs(x(1:end-WL_delay)), 20*log10(abs(transfer_no_PD)), 'o');
hold on;
plot(abs(x(1:end-WL_delay)), 20*log10(abs(transfer_with_MP_PD)), '*');
legend('Without DPD','with MP PD')
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
pspectrum(y_d(mem_depth+1:end), Fs)
hold on
pspectrum(y(mem_depth+1:end), Fs)
hold on
pspectrum(y_MP(mem_depth+1:end), Fs)
legend('y_d','no PD','with MP PD')
