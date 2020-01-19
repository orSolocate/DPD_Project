clc; clear; close all;                                     %start freash

%% General User inputs
%amp_data          = 'simrfV2_powamp_dpd_data.mat';
%signal            = '100MHzLTE.mat';
signal            = 'testsignal.mat';
in_sig            = 'inDataPA.mat';
out_sig           = 'outDataPA.mat';
model_run_period  = 50000 ;                                %num of samples to run in the model
start_pos_sig     = 50001;
end_pos_sig       = start_pos_sig+model_run_period-1;

%% Loads
load(in_sig);
load(out_sig);
inDataPA  = x';
outDataPA = y';
load(signal);                                              %load signal
%load(amp_data);                                            %load amp data

%% Model User inputs
iterations_num_MP  = 2;
resistor           = 50;
miu_MP             = -0.05;

mem_depth = 2 ;                                            %M in the MP model
mem_deg   = 5 ;                                            %K in the MP model

%% Initialize and Model
avg_gain         = abs(mean(outDataPA./inDataPA));         
%x                = waveform(start_pos_sig:end_pos_sig);
x                = x(start_pos_sig:end_pos_sig)./norm(x,2);
y_d              = x.*avg_gain;
AMP_coef_Matrix  = Get_coef_MP(inDataPA, outDataPA, mem_deg, mem_depth);
PD_coef_Matrix   = Get_coef_MP(outDataPA, inDataPA, mem_deg, mem_depth);

%dir="Experiments/MP_"+iterations_num_MP+"_iteration_"+mem_deg+"_deg_"+mem_depth+"_depth_"+miu_MP+"_miu"
%mkdir(dir)

%% initial calculations
x_opt_MP  = [zeros(mem_depth,1); PD_MP(y_d, PD_coef_Matrix, mem_deg, mem_depth)];

y_MP     = [zeros(mem_depth,1); Get_model_output_MP(AMP_coef_Matrix, x_opt_MP, mem_deg, mem_depth)];
y_in_MP  = [zeros(mem_depth,1); Get_model_output_MP(AMP_coef_Matrix, x, mem_deg, mem_depth)];

%phase correction
%x_opt_MP  = ifft(fft(x_opt_MP).*exp(-phdiffmeasure(y_MP, x_opt_MP)*1i));

%% Sub calculations
transfer_no_PD        = y_in_MP./x;
transfer_with_MP_PD   = y_MP./x_opt_MP;

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
    PD_coef_Matrix = Update_coef_MP(y_d, x_opt_MP, PD_coef_Matrix, mem_deg, mem_depth, miu_MP);
    
    %get the updated optimal input
    x_opt_MP       = [zeros(mem_depth,1); PD_MP(y_d, PD_coef_Matrix, mem_deg, mem_depth)];
    
    %phase correction
    %x_opt_MP       = ifft(fft(x_opt_MP).*exp(-phdiffmeasure(y_MP, x_opt_MP)*1i));
    
    %get the Amp output
    y_MP  = [zeros(mem_depth,1); Get_model_output_MP(AMP_coef_Matrix, x_opt_MP, mem_deg, mem_depth)];
    
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
    
    if(0)
    if((error_per_iter_MP(ll) == error_per_iter_MP(ll-1)))
        miu_MP = 0.5*miu_MP; 
    end
    if((error_per_iter_MP(ll) > error_per_iter_MP(ll-1)))
        miu_MP = -0.5*miu_MP; 
    end
    if((error_per_iter_MP(ll) < error_per_iter_MP(ll-1)))
        if(miu_MP < -0.001)
            miu_MP = 1.1*miu_MP; 
        end
    end
    end    
end

%% Sub calculations
transfer_no_PD        = y_in_MP./x;
transfer_with_MP_PD   = y_MP./x_opt_MP;

%% Ploting final results
figure;
plot(abs(x), 20*log10(abs(transfer_no_PD)), 'o');
hold on;
plot(abs(x_opt_MP), 20*log10(abs(transfer_with_MP_PD)), '*');
legend('Without DPD in MP','with MP PD')
xlabel('abs(input)');
ylabel('Gain');
title('Linearity Comparison');
%savefig(dir+"/linearity.fig")

figure
plot(error_per_iter_MP(2:end))
legend('error per iteration MP')
xlabel('iteration #');
ylabel('error on y (norm-2)');
title("Error (iteration)")
%savefig(dir+"/error_y.fig")

figure
plot(phase_per_iter_MP)
legend('phase per iteration MP')
xlabel('iteration #');
ylabel('phase difference [degrees]');
title("Phase (iteration)")
%savefig(dir+"/phase.fig")

figure
plot(avg_power_per_iter_MP)
legend('avg power per iteration MP')
xlabel('iteration #');
ylabel('average input power [dBm]');
title("Average Power (iteration)")
%savefig(dir+"/in_power.fig")

figure
plot(RMS_in_per_iter_MP)
legend('RMS in per iteration MP')
xlabel('iteration #');
ylabel('RMS in [dBm]');
title("RMS in(iteration)")
%savefig(dir+"/rms_in.fig")
%%
%[y_MP, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_1(x_opt_MP);
%[y_in_MP, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_1(x);
figure
pspectrum(y_d(mem_depth+1:end), 200e6)
hold on
pspectrum(y_in_MP(mem_depth+1:end), 200e6)
hold on
pspectrum(y_MP(mem_depth+1:end), 200e6)
legend('optimal output','no PD','with MP PD')
%title("Signal Y After "+iterations_num_MP+" iterations")
%savefig(dir+"/y.fig")