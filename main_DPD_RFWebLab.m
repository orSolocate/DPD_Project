clc; clear; close all;                                     %start freash

%% General User inputs
signal            = '100MHzLTE.mat';
model_run_period  = 100000;                                %num of samples to run in the model
start_pos_sig     = 250001;
end_pos_sig       = start_pos_sig+model_run_period-1;
coef_update_time  = 10000000;
phase_ud_time     = 10000000;
phase_diff        = 0.1;


%% RF-WEBLAB User inputs
memory_deg_weblab   = 3;                                   %K in the MP model
memory_depth_weblab = 5;                                   %M in the MP model
Ma_WL = 9; Mb_WL = 2; Mc_WL = 0;                           %memory deg for GMP model
Ka_WL = 8; Kb_WL = 5; Kc_WL = 0;                           %non-linearity deg for GMP
P_WL = 5; Q_WL = 0;                                        %cross-terms for GMP


%% RF-WEBLAB Initialize
load(signal);                                              %load signal
x                = waveform(start_pos_sig:end_pos_sig);
[y, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(x);
weblab_delay           = finddelay(x, y);
orders_WL              = [Ma_WL, Ka_WL, Mb_WL, Kb_WL, P_WL, Mc_WL, Kc_WL, Q_WL];

if(weblab_delay >=0)
    avg_gain_weblab    = abs(mean(y(weblab_delay+1:end)./x(1:end-weblab_delay)));
    PD_coef_Matrix_WL  = Get_coef_MP(y(weblab_delay+1:end)', x(1:end-weblab_delay)', memory_deg_weblab, memory_depth_weblab);
    PD_coef_Vec_GMP_WL = Get_coef_GMP(y(weblab_delay+1:end)', x(1:end-weblab_delay)', orders_WL);
end
if(weblab_delay < 0)
    avg_gain_weblab    = abs(mean(y(1:end-weblab_delay)./x(weblab_delay+1:end)));
    PD_coef_Matrix_WL  = Get_coef_MP(y(1:end-weblab_delay)', x(weblab_delay+1:end)', memory_deg_weblab, memory_depth_weblab);
    PD_coef_Vec_GMP_WL = Get_coef_GMP(y(1:end-weblab_delay)', x(weblab_delay+1:end)', orders_WL);
end
y_d_WL                 = x.*avg_gain_weblab;
x_opt_weblab_MP        = zeros(model_run_period,1);
x_opt_weblab_GMP       = zeros(model_run_period,1);
first_n_GMP_WL         = 1+max([orders_WL(1), orders_WL(3)+orders_WL(5), orders_WL(6)-1]);


%% RF-WEBLAB MP Run
n = memory_depth_weblab+1;
while (n<length(y_d_WL))                                         %run the system as live
    x_opt_weblab_MP(n) = PD_MP(y_d_WL(n-memory_depth_weblab:n), PD_coef_Matrix_WL, memory_deg_weblab, memory_depth_weblab);
    n                  = n+1;
end

%% RF-WEBLAB GMP Run
n = first_n_GMP_WL;
while (n<length(y_d_WL)-Q_WL)                                    %run the system as live
    x_opt_weblab_GMP(n) = PD_GMP(y_d_WL(n-first_n_GMP_WL+1:n+Q_WL), PD_coef_Vec_GMP_WL, orders_WL);
    n                   = n+1;
end


%% RF-WEBLAB Ploting
[y_RFW_MP, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(x_opt_weblab_MP);
[y_RFW_GMP, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_1(x_opt_weblab_GMP);
transfer_no_PD_weblab         = y./x;
transfer_with_MP_PD_weblab    = y_RFW_MP./x_opt_weblab_MP;
transfer_with_GMP_PD_weblab   = y_RFW_GMP./x_opt_weblab_GMP;

figure
pspectrum(y_d_WL(memory_depth_weblab+1:end))
hold on
pspectrum(y)
hold on
pspectrum(y_RFW_MP(memory_depth_weblab+1:end))
hold on
pspectrum(y_RFW_GMP)
legend('y_d','no PD','with MP PD','with GMP PD')

figure;
plot(abs(y_d_WL./avg_gain_weblab), 20*log10(abs(transfer_no_PD_weblab)), 'o');
hold on;
plot(abs(y_d_WL./avg_gain_weblab), 20*log10(abs(transfer_with_MP_PD_weblab)), '+');
hold on;
plot(abs(y_d_WL./avg_gain_weblab), 20*log10(abs(transfer_with_GMP_PD_weblab)), '*');
legend('Without DPD','with MP PD','with GMP PD')
xlabel('abs(input)');
ylabel('Magnitude');
title('Linearity Comparison');