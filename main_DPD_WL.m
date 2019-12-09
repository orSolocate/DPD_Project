clc; clear; close all;                                     %start freash

%% General User inputs
signal            = 'testsignal.mat';
model_run_period  = 50000 ;                                %num of samples to run in the model
start_pos_sig     = 25001;
end_pos_sig       = start_pos_sig+model_run_period-1;
iterations_num    = 20;
miu_MP            = 0.03;
miu_GMP           = 0.03;

%% RF-WEBLAB User inputs
mem_deg   = 3;                                              %K in the MP model
mem_depth = 5;                                              %M in the MP model
Ma_WL = 9 ; Mb_WL = 2 ; Mc_WL = 2 ;                         %memory deg for GMP model
Ka_WL = 8 ; Kb_WL = 5 ; Kc_WL = 5 ;                         %non-linearity deg for GMP
             P_WL = 10;  Q_WL = 10;                         %cross-terms for GMP

%% RF-WEBLAB Initialize and Model
load(signal);                                              %load signal
x                      = x(start_pos_sig:end_pos_sig)./(norm(x,2));
[y, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(x);
WL_delay               = finddelay(x, y);
orders_WL              = [Ma_WL, Ka_WL, Mb_WL, Kb_WL, P_WL, Mc_WL, Kc_WL, Q_WL];

if(WL_delay >=0)
    avg_gain           = abs(mean(y(WL_delay+1:end)./x(1:end-WL_delay)));
    PD_coef_Matrix_WL  = Get_coef_MP(y(WL_delay+1:end)', x(1:end-WL_delay)', mem_deg, mem_depth);
    PD_coef_Vec_GMP_WL = Get_coef_GMP(y(WL_delay+1:end)', x(1:end-WL_delay)', orders_WL);
end
if(WL_delay < 0)
    avg_gain           = abs(mean(y(1:end-WL_delay)./x(WL_delay+1:end)));
    PD_coef_Matrix_WL  = Get_coef_MP(y(1:end-WL_delay)', x(WL_delay+1:end)', mem_deg, mem_depth);
    PD_coef_Vec_GMP_WL = Get_coef_GMP(y(1:end-WL_delay)', x(WL_delay+1:end)', orders_WL);
end
y_d_WL                 = x.*avg_gain;
first_n_GMP_WL         = 1+max([orders_WL(1), orders_WL(3)+orders_WL(5), orders_WL(6)-1]);

%% RF-WEBLAB initial calculations
x_opt_weblab_MP  = [zeros(mem_depth,1); PD_MP(y_d_WL, PD_coef_Matrix_WL, mem_deg, mem_depth)];
x_opt_weblab_GMP = [zeros(first_n_GMP_WL-1,1); PD_GMP(y_d_WL, PD_coef_Vec_GMP_WL,orders_WL); zeros(orders_WL(8),1)];
[y_RFW_MP, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(x_opt_weblab_MP);
[y_RFW_GMP, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_1(x_opt_weblab_GMP);
transfer_no_PD_weblab         = y./x;
transfer_with_MP_PD_weblab    = y_RFW_MP./x_opt_weblab_MP;
transfer_with_GMP_PD_weblab   = y_RFW_GMP./x_opt_weblab_GMP;

%% RF-WEBLAB Ploting
figure
pspectrum(y_d_WL(mem_depth+1:end))
hold on
pspectrum(y)
hold on
pspectrum(y_RFW_MP(mem_depth+1:end))
hold on
pspectrum(y_RFW_GMP)
legend('y_d','no PD','with MP PD','with GMP PD')

figure;
plot(abs(x), 20*log10(abs(transfer_no_PD_weblab)), 'o');
hold on;
plot(abs(x_opt_weblab_MP), 20*log10(abs(transfer_with_MP_PD_weblab)), '+');
hold on;
plot(abs(x_opt_weblab_GMP), 20*log10(abs(transfer_with_GMP_PD_weblab)), '*');
legend('Without DPD','with MP PD','with GMP PD')
xlabel('abs(input)');
ylabel('Magnitude');
title('Linearity Comparison');

%% Iterate
error_per_iter_MP            = ones(iterations_num+1,1).*inf;
error_per_iter_GMP           = ones(iterations_num+1,1).*inf;
ll                           = 1;
while (ll<=iterations_num)
    
    x_opt_tilda_weblab_MP    = [zeros(mem_depth,1); PD_MP(y_RFW_MP, PD_coef_Matrix_WL, mem_deg, mem_depth)];
    x_opt_tilda_weblab_GMP   = [zeros(first_n_GMP_WL-1,1); PD_GMP(y_RFW_GMP, PD_coef_Vec_GMP_WL,orders_WL); zeros(orders_WL(8),1)];
    error_x_opts_MP          = x_opt_weblab_MP - x_opt_tilda_weblab_MP;
    error_x_opts_GMP         = x_opt_weblab_GMP - x_opt_tilda_weblab_GMP;
    error_per_iter_MP(ll+1)  = norm(error_x_opts_MP,2);
    error_per_iter_GMP(ll+1) = norm(error_x_opts_GMP,2);
    %err_MP                   = (y_RFW_MP./avg_gain) - x;
    %err_MP                   = (y_RFW_GMP./avg_gain) - x;
    
    %get the updated coefficients for the model
    PD_coef_Matrix_WL  = Update_coef_MP(error_x_opts_MP, y_d_WL, PD_coef_Matrix_WL, mem_deg, mem_depth, miu_MP);
    PD_coef_Vec_GMP_WL = Update_coef_GMP(error_x_opts_GMP, y_d_WL, PD_coef_Vec_GMP_WL, orders_WL, miu_GMP);
    
    %get the updated optimal input
    x_opt_weblab_MP        = [zeros(mem_depth,1); PD_MP(y_d_WL, PD_coef_Matrix_WL, mem_deg, mem_depth)];
    x_opt_weblab_GMP       = [zeros(first_n_GMP_WL-1,1); PD_GMP(y_d_WL, PD_coef_Vec_GMP_WL,orders_WL); zeros(orders_WL(8),1)];
    
    %phase correction
    x_opt_weblab_MP        = ifft(fft(x_opt_weblab_MP).*exp(-phdiffmeasure(y_RFW_MP, x_opt_weblab_MP)*1i));
    x_opt_weblab_GMP       = ifft(fft(x_opt_weblab_GMP).*exp(-phdiffmeasure(y_RFW_GMP, x_opt_weblab_GMP)*1i));
    
    %get the Amp output
    [y_RFW_MP, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(x_opt_weblab_MP);
    [y_RFW_GMP, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_1(x_opt_weblab_GMP);
    
    %spectrum plot
    figure
    pspectrum(y_d_WL(mem_depth+1:end))
    hold on
    pspectrum(y)
    hold on
    pspectrum(y_RFW_MP(mem_depth+1:end))
    hold on
    pspectrum(y_RFW_GMP)
    legend('y_d','no PD','with MP PD','with GMP PD')  
    
    disp(['iteration num: ',num2str(ll), ' is done.'])
    ll = ll+1;
    
    if(error_per_iter_MP(ll)>error_per_iter_MP(ll-1))
        miu_MP = -miu_MP;
    end
    
    if(error_per_iter_GMP(ll)>error_per_iter_GMP(ll-1))
        miu_GMP = -miu_GMP;
    end
    
end

%%
figure
plot(error_per_iter_MP(2:end))
hold on
plot(error_per_iter_GMP(2:end))
legend('error per iteration MP','error per iteration GMP')

