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
iteration = 3                                               

%% RF-WEBLAB Initialize
load(signal);                                              %load signal
x                = waveform(start_pos_sig:end_pos_sig);
[y, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(x);

%% ILC_Scheme
orders=[memory_deg_weblab,memory_depth_weblab];
[x_opt, error_vec_plot] = ILC_Scheme_RFWebLab(x,y, orders,100);
figure; 
plot(1:length(error_vec_plot), error_vec_plot);
xlabel('Iterations');
ylabel('Error');
title('Error vs. Iteration Number');
grid;
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','GridAlpha'),'GridAlpha',1);




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