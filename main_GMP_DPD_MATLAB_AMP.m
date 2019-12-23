clc; clear; close all;                                     %start freash

%% General User inputs
amp_data          = 'simrfV2_powamp_dpd_data.mat';
signal            = '100MHzLTE.mat';
model_run_period  = 60000 ;                                %num of samples to run in the model
start_pos_sig     = 1;
end_pos_sig       = start_pos_sig+model_run_period-1;

%in comment below: not used for now - maybe will be useful later
%coef_update_time  = 10000000;
%phase_ud_time     = 10000000;
%phase_diff        = 0.1;

%% Loads
load(signal);                                              %load signal
load(amp_data);                                            %load amp data
%% Model User inputs
iterations_num_GMP = 20;
resistor           = 50;
miu_GMP            = -0.009;

Ma = 2 ; Mb = 1 ; Mc = 1 ;                                 %memory depth for GMP model
Ka = 5 ; Kb = 4 ; Kc = 4 ;                                 %non-linearity deg for GMP
          P = 6 ;  Q = 6 ;                                 %cross-terms for GMP

%% Initialize and Model
orders           = [Ma, Ka, Mb, Kb, P, Mc, Kc, Q];
first_n_GMP      = 1+max([orders(1), orders(3)+orders(5), orders(6)-1]);
avg_gain         = abs(mean(outDataPA./inDataPA)); %in matlab AMP there is no delay        
x                = waveform(start_pos_sig:end_pos_sig);
y_d              = x.*avg_gain;
AMP_coef_GMP_vec = Get_coef_GMP(inDataPA, outDataPA, orders);
PD_coef_GMP_vec  = Get_coef_GMP(outDataPA, inDataPA, orders);

%% initial calculations
x_opt_GMP = [zeros(first_n_GMP-1,1); PD_GMP(y_d, PD_coef_GMP_vec,orders); zeros(orders(8),1)];

y_GMP    = [zeros(first_n_GMP-1,1); Get_model_output_GMP(AMP_coef_GMP_vec, x_opt_GMP, orders); zeros(orders(8),1)];
y_in_GMP = [zeros(first_n_GMP-1,1); Get_model_output_GMP(AMP_coef_GMP_vec, x, orders); zeros(orders(8),1)];

%phase correction
%x_opt_GMP = ifft(fft(x_opt_GMP).*exp(-phdiffmeasure(y_GMP, x_opt_GMP)*1i));

%% Sub calculations
transfer_no_PD_in_GMP = y_in_GMP./x;
transfer_with_GMP_PD  = y_GMP./x_opt_GMP;

%% Ploting initial state
figure
pspectrum(y_d(orders(1)+1:end))
hold on
pspectrum(y_in_GMP(orders(1)+1:end))
hold on
pspectrum(y_GMP(orders(1)+1:end))
legend('y_d','no PD in GMP','with GMP PD')

figure;
plot(abs(x), 20*log10(abs(transfer_no_PD_in_GMP)), 'o');
hold on;
plot(abs(x_opt_GMP), 20*log10(abs(transfer_with_GMP_PD)), '*');
legend('Without DPD in GMP','with GMP PD')
xlabel('abs(input)');
ylabel('Gain');
title('Linearity Comparison');
%% Iterate pre
error_per_iter_GMP        = ones(iterations_num_GMP+1,1).*inf;

phase_per_iter_GMP        = zeros(iterations_num_GMP+1,1);
phase_per_iter_GMP(1)     = phdiffmeasure(y_GMP, x_opt_GMP)*180/pi;

avg_power_per_iter_GMP    = zeros(iterations_num_GMP+1,1);
avg_power_per_iter_GMP(1) = 10*log10( norm(x_opt_GMP)^2/resistor/length(x)) + 30;

PAPR_per_iter_GMP         = zeros(iterations_num_GMP+1,1);
PAPR_per_iter_GMP(1)      = 20*log10(max(abs(x_opt_GMP))*sqrt(length(x_opt_GMP))/norm(x_opt_GMP));

RMS_in_per_iter_GMP       = zeros(iterations_num_GMP+1,1);
RMS_in_per_iter_GMP(1)    = avg_power_per_iter_GMP(1);

%% Iterate GMP
ll = 1;
while (ll<=iterations_num_GMP)
    
    y_GMP                    = ifft(fft(y_GMP).*exp(-phdiffmeasure(y_d, y_GMP)*1i));
    output_err_GMP           = y_d - y_GMP;
    error_per_iter_GMP(ll+1) = norm(output_err_GMP,2);
    
    %get the updated coefficients for the model
    PD_coef_GMP_vec = Update_coef_GMP(y_d, x_opt_GMP, PD_coef_GMP_vec, orders, miu_GMP);
    
    %get the updated optimal input
    x_opt_GMP       = [zeros(first_n_GMP-1,1); PD_GMP(y_d, PD_coef_GMP_vec,orders); zeros(orders(8),1)];

    %phase correction
    %x_opt_GMP       = ifft(fft(x_opt_GMP).*exp(-phdiffmeasure(y_GMP, x_opt_GMP)*1i));

    %get the Amp output
    y_GMP = [zeros(first_n_GMP-1,1); Get_model_output_GMP(AMP_coef_GMP_vec, x_opt_GMP, orders); zeros(orders(8),1)];

    %spectrum plot every 2nd iteration
    if(0 == mod(ll,2))
        figure
        pspectrum(y_d(orders(1)+1:end))
        hold on
        pspectrum(y_in_GMP(orders(1)+1:end))
        hold on
        pspectrum(y_GMP(orders(1)+1:end))
        legend('y_d','no PD','with GMP PD')
    end
       
    disp(['GMP iteration num: ',num2str(ll), ' is done.'])
    ll = ll+1;
    
    phase_per_iter_GMP(ll)     = phdiffmeasure(y_GMP, x_opt_GMP)*180/pi;
    avg_power_per_iter_GMP(ll) = 10*log10( norm(x_opt_GMP)^2/resistor/length(x)) + 30;
    PAPR_per_iter_GMP(ll)      = 20*log10(max(abs(x_opt_GMP))*sqrt(length(x_opt_GMP))/norm(x_opt_GMP));
    RMS_in_per_iter_GMP(ll)    = RMS_in_per_iter_GMP(ll-1) - (PAPR_per_iter_GMP(ll) - PAPR_per_iter_GMP(ll-1));
    
    if((error_per_iter_GMP(ll) == error_per_iter_GMP(ll-1)))
        miu_GMP = 0.5*miu_GMP; 
    end
    if((error_per_iter_GMP(ll) > error_per_iter_GMP(ll-1)))
        miu_GMP = -0.5*miu_GMP; 
    end
    if((error_per_iter_GMP(ll) < error_per_iter_GMP(ll-1)))
        if(miu_GMP < -0.001)
            miu_GMP = 1.1*miu_GMP; 
        end
    end
    
end

%% Sub calculations
transfer_no_PD_in_GMP = y_in_GMP./x;
transfer_with_GMP_PD  = y_GMP./x_opt_GMP;

%% Ploting final results
figure;
plot(abs(x), 20*log10(abs(transfer_no_PD_in_GMP)), 'o');
hold on;
plot(abs(x_opt_GMP), 20*log10(abs(transfer_with_GMP_PD)), '*');
legend('Without DPD in GMP','with GMP PD')
xlabel('abs(input)');
ylabel('Gain');
title('Linearity Comparison');

figure
plot(error_per_iter_GMP(2:end))
legend('error per iteration GMP')
xlabel('iteration #');
ylabel('error on y (norm-2)');

figure
plot(phase_per_iter_GMP)
legend('phase per iteration GMP')
xlabel('iteration #');
ylabel('phase diff (deg)');

figure
plot(avg_power_per_iter_GMP)
legend('avg power per iteration GMP')
xlabel('iteration #');
ylabel('average input power (dBm)');

figure
plot(RMS_in_per_iter_GMP)
legend('RMS in per iteration GMP')
xlabel('iteration #');
ylabel('RMS in');

figure
pspectrum(y_d(orders(1)+1:end))
hold on
pspectrum(y_in_GMP(orders(1)+1:end))
hold on
pspectrum(y_GMP(orders(1)+1:end))
legend('y_d','no PD in GMP','with GMP PD')
