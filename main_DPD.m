clc; clear; close all;                                     %start freash

%% General User inputs
amp_data          = 'simrfV2_powamp_dpd_data.mat';
signal            = '100MHzLTE.mat';
model_run_period  = 100000;                                %num of samples to run in the model
start_pos_sig     = 250001;
end_pos_sig       = start_pos_sig+model_run_period-1;
coef_update_time  = 10000000;
phase_ud_time     = 10000000;
phase_diff        = 0.1;

%% MATLAB AMP User inputs
memory_deg        = 5;                                     %K in the MP model
memory_depth      = 2;                                     %M in the MP model
Ma = 2; Mb = 1; Mc = 1;                                    %memory deg for GMP model
Ka = 5; Kb = 7; Kc = 7;                                    %non-linearity deg for GMP
P = 5; Q = 5;                                              %cross-terms for GMP

%% RF-WEBLAB User inputs
memory_deg_weblab   = 3;                                   %K in the MP model
memory_depth_weblab = 5;                                   %M in the MP model
Ma_WL = 7; Mb_WL = 2; Mc_WL = 2;                           %memory deg for GMP model
Ka_WL = 8; Kb_WL = 5; Kc_WL = 5;                           %non-linearity deg for GMP
P_WL = 12; Q_WL = 12;                                        %cross-terms for GMP

%% Initialize
load(amp_data);                                            %load amp data
load(signal);                                              %load signal
avg_gain         = abs(mean(outDataPA./inDataPA));         %set the avg gain
x_opt_MP         = zeros(model_run_period,1);              %optimal input to the amp
x_opt_GMP        = zeros(model_run_period,1);              %optimal input to the amp 
output_MP        = zeros(model_run_period,1);              %set an output vector
output_GMP       = zeros(model_run_period,1);              %set an output vector
x                = waveform(start_pos_sig:end_pos_sig);
y_d              = x.*avg_gain;                            %get the desired output for model run
orders           = [Ma, Ka, Mb, Kb, P, Mc, Kc, Q];
AMP_coef_Matrix  = Get_coef_MP(inDataPA, outDataPA, memory_deg, memory_depth);
PD_coef_Matrix   = Get_coef_MP(outDataPA, inDataPA, memory_deg, memory_depth);
AMP_coef_GMP_vec = Get_coef_GMP(inDataPA, outDataPA, orders);
PD_coef_GMP_vec  = Get_coef_GMP(outDataPA, inDataPA, orders);
first_n_GMP      = 1+max([orders(1), orders(3)+orders(5), orders(6)-1]);

%% RF-WEBLAB Initialize
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

%% MATLAB AMP MP Run
n = memory_depth+1;
while (n<length(y_d))                                      %run the system as live
    %y_for_PD  = Get_y_for_PD_MP(y_d, output, n, memory_depth);
    x_opt_MP(n)  = PD_MP(y_d(n-memory_depth:n), PD_coef_Matrix, memory_deg, memory_depth);
    %if (0 == mod(n,phase_ud_time))
    %    if(abs(phdiffmeasure(output, x_opt)) >= phase_diff)
    %    x_opt     = ifft(fft(x_opt).*exp(-phdiffmeasure(output, x_opt)*1i));
    %    end
    %end
    output_MP(n) = Get_model_output_MP(AMP_coef_Matrix, x_opt_MP(n-memory_depth:n), memory_deg, memory_depth);
    %if (0 == mod(n,coef_update_time))
    %    PD_coef_Matrix  = Get_coef_MP([output(1:n)' outDataPA], [x_opt(1:n)' inDataPA], memory_deg, memory_depth);
    %end
    n = n+1;
end

%% MATLAB AMP GMP Run
n = first_n_GMP;
while (n<length(y_d)-Q)                                    %run the system as live
    x_opt_GMP(n)  = PD_GMP(y_d(n-first_n_GMP+1:n+Q), PD_coef_GMP_vec, orders);
    output_GMP(n) = Get_model_output_GMP(AMP_coef_GMP_vec, x_opt_GMP(n-first_n_GMP+1:n+Q), orders);
    n = n+1;
end

%% MATLAB AMP Ploting 
y_no_PD = Get_model_output_MP(AMP_coef_Matrix, y_d./avg_gain, memory_deg, memory_depth);
transfer_no_PD   = y_no_PD./(y_d(1:length(y_no_PD))./avg_gain);
transfer_with_PD = output_MP(memory_depth+1:end)./(y_d(1:length(y_no_PD))./avg_gain);

figure
pspectrum(y_d(memory_depth+1:end))
hold on
pspectrum(y_no_PD)
hold on
pspectrum(output_MP(memory_depth+1:end))
hold on
pspectrum(output_GMP)
legend('y_d','no PD','with MP PD','with GMP PD')

figure;
plot(abs(y_d(1:1:length(y_no_PD))./avg_gain), 20*log10(abs(transfer_no_PD(1:1:end))), 'o');
hold on;
plot(abs(y_d(1:1:length(y_no_PD))./avg_gain), 20*log10(abs(transfer_with_PD(1:1:end))), '+');
legend('Without DPD','With DPD')
xlabel('abs(input)');
ylabel('Magnitude');
title('Linearity Comparison');

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
transfer_with_MP_PD_weblab    = y_RFW_MP./x;
transfer_with_GMP_PD_weblab   = y_RFW_GMP./x;

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
plot(abs(y_d_WL./avg_gain), 20*log10(abs(transfer_no_PD_weblab)), 'o');
hold on;
plot(abs(y_d_WL./avg_gain), 20*log10(abs(transfer_with_MP_PD_weblab)), '+');
hold on;
plot(abs(y_d_WL./avg_gain), 20*log10(abs(transfer_with_GMP_PD_weblab)), '*');
legend('Without DPD','with MP PD','with GMP PD')
xlabel('abs(input)');
ylabel('Magnitude');
title('Linearity Comparison');