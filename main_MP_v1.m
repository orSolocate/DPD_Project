clc; clear; close all;                                     %start freash

%% General User inputs
%signal            = '100MHzLTE.mat';
signal            = 'testsignal.mat';
model_run_period  = 50000 ;                                %num of samples to run in the model
start_pos_sig     = 25001;
end_pos_sig       = start_pos_sig+model_run_period-1;
iterations_num    = 2;
miu_MP            = 0.05;

%% RF-WEBLAB User inputs
mem_deg   = 3 ;                                             %K in the MP model
mem_depth = 5 ;                                             %M in the MP model

%% RF-WEBLAB Initialize and Model
load(signal);                                              %load signal
%x                      = waveform(start_pos_sig:end_pos_sig);
x                      = x(start_pos_sig:end_pos_sig)./(norm(x,2));
[y, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(x);
WL_delay               = finddelay(x, y);

if(WL_delay >=0)
    avg_gain           = abs(mean(y(WL_delay+1:end)./x(1:end-WL_delay)));
    PD_coef_Matrix_WL  = Get_coef_MP(y(WL_delay+1:end)', x(1:end-WL_delay)', mem_deg, mem_depth);
end
if(WL_delay < 0)
    avg_gain           = abs(mean(y(1:end-WL_delay)./x(WL_delay+1:end)));
    PD_coef_Matrix_WL  = Get_coef_MP(y(1:end-WL_delay)', x(WL_delay+1:end)', mem_deg, mem_depth);
end
y_d_WL                 = x.*avg_gain;

%% RF-WEBLAB initial calculations
x_opt_weblab_MP  = [zeros(mem_depth,1); PD_MP(y_d_WL, PD_coef_Matrix_WL, mem_deg, mem_depth)];
[y_RFW_MP, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(x_opt_weblab_MP);
transfer_no_PD_weblab         = y./x;
transfer_with_MP_PD_weblab    = y_RFW_MP./x_opt_weblab_MP;

%% RF-WEBLAB Ploting
figure
pspectrum(y_d_WL(mem_depth+1:end))
hold on
pspectrum(y)
hold on
pspectrum(y_RFW_MP(mem_depth+1:end))
legend('y_d','no PD','with MP PD')

figure;
plot(abs(x), 20*log10(abs(transfer_no_PD_weblab)), 'o');
hold on;
plot(abs(x_opt_weblab_MP), 20*log10(abs(transfer_with_MP_PD_weblab)), '+');
legend('Without DPD','with MP PD')
xlabel('abs(input)');
ylabel('Magnitude');
title('Linearity Comparison');

%% Iterate
error_per_iter_MP            = ones(iterations_num+1,1).*inf;
ll                           = 1;
while (ll<=iterations_num)
    
    x_opt_tilda_weblab_MP    = [zeros(mem_depth,1); PD_MP(y_RFW_MP, PD_coef_Matrix_WL, mem_deg, mem_depth)];
    error_x_opts_MP          = x_opt_weblab_MP - x_opt_tilda_weblab_MP;
    error_per_iter_MP(ll+1)  = norm(error_x_opts_MP,2);
    %err_MP                   = (y_RFW_MP./avg_gain) - x;
    
    %get the updated coefficients for the model
    PD_coef_Matrix_WL  = Update_coef_MP(error_x_opts_MP, y_d_WL, PD_coef_Matrix_WL, mem_deg, mem_depth, miu_MP);
    
    %get the updated optimal input
    x_opt_weblab_MP        = [zeros(mem_depth,1); PD_MP(y_d_WL, PD_coef_Matrix_WL, mem_deg, mem_depth)];
    
    %phase correction
    x_opt_weblab_MP        = ifft(fft(x_opt_weblab_MP).*exp(-phdiffmeasure(y_RFW_MP, x_opt_weblab_MP)*1i));
    
    %get the Amp output
    [y_RFW_MP, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(x_opt_weblab_MP);
    
    %spectrum plot
    figure
    pspectrum(y_d_WL(mem_depth+1:end))
    hold on
    pspectrum(y)
    hold on
    pspectrum(y_RFW_MP(mem_depth+1:end))
    legend('y_d','no PD','with MP PD')  
    
    disp(['iteration num: ',num2str(ll), ' is done.'])
    ll = ll+1;
    
    %if(error_per_iter_MP(ll)>error_per_iter_MP(ll-1))
    %    miu_MP = -miu_MP;
    %end

    
end

%%
figure
plot(error_per_iter_MP(2:end))
legend('error per iteration MP')

