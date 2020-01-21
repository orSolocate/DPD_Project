clc; clear; close all; 

%% General User inputs
amp_data          = 'simrfV2_powamp_dpd_data.mat';
signal            = 'testsignal.mat';
model_run_period  = 50000 ;                                %num of samples to run in the model
start_pos_sig     = 1;
end_pos_sig       = start_pos_sig+model_run_period-1;

%% Model User inputs
miu_GMP = 0.005;
Ma = 5 ; Mb = 5 ; Mc = 5 ;                                 %memory depth for GMP model
Ka = 9 ; Kb = 5 ; Kc = 5 ;                                 %non-linearity deg for GMP
          P = 3 ;  Q = 3 ;                                 %cross-terms for GMP

orders      = [Ma, Ka, Mb, Kb, P, Mc, Kc, Q];
first_n_GMP = 1+max([orders(1), orders(3)+orders(5), orders(6)-1]);

%% Loads
load(signal);                                              %load signal
load(amp_data);                                            %load amp data
x            = x(start_pos_sig:end_pos_sig)./(2*norm(x,2));
AMP_coef_Vec = Get_coef_GMP(inDataPA, outDataPA, orders);
G            = mean(abs(outDataPA))/mean(abs(inDataPA));%abs(mean(y./x));
y            = [zeros(first_n_GMP-1,1); Get_model_output_GMP(AMP_coef_Vec, x.*G, orders); zeros(orders(8),1)];

%% Initialize and Model
y    = ifft(fft(y).*exp(-phdiffmeasure(x, y)*1i));
coef = Get_coef_GMP(y', x', orders);

%% Itarations pre
k                = 1;
error_decreases  = true;
DISPLAY_FACTOR   = 10;
max_iter         = 10;
error_lim        = 1e-5;
iter_error       = ones(1,1).*inf;

%% First iteration
z             = [zeros(first_n_GMP-1,1);PD_GMP(x, coef, orders); zeros(orders(8),1)];
z_hat         = [zeros(first_n_GMP-1,1);PD_GMP(y, coef, orders); zeros(orders(8),1)];
error         = z - z_hat;
iter_error(k) = norm(error,2);

%% Itarate
while((iter_error(k) > error_lim) && (error_decreases) && (k < max_iter))
    k = k + 1;
    
    delay    = finddelay(error,y);
    coef_hat = Get_coef_GMP(y', error', orders);
    coef     = (coef + miu_GMP.*coef_hat);  
    
    z               = [zeros(first_n_GMP-1,1);PD_GMP(x, coef, orders); zeros(orders(8),1)];
    z_hat           = [zeros(first_n_GMP-1,1);PD_GMP(y, coef, orders); zeros(orders(8),1)];
    error           = z - z_hat;
    iter_error(k)   = norm(error,2);
    error_decreases = iter_error(k) <= iter_error(k-1);
    if(rem(k,DISPLAY_FACTOR) == 0)
        disp(['Iteration #: ',num2str(k), ' - Error is: ',num2str(iter_error(k))]);
    end
 
end

%% Plots
figure; plot(iter_error);

z_new    = [zeros(first_n_GMP-1,1);PD_GMP(x.*G, coef, orders); zeros(orders(8),1)];

y_GMP    = [zeros(first_n_GMP-1,1); Get_model_output_GMP(AMP_coef_Vec, z_new.*G, orders); zeros(orders(8),1)].*G;
y_no_GMP = [zeros(first_n_GMP-1,1); Get_model_output_GMP(AMP_coef_Vec, x.*G, orders); zeros(orders(8),1)];

figure; pspectrum(y_GMP, 200e6); hold on; pspectrum(y_no_GMP, 200e6); legend('y_{GMP}','y_{no GMP}');
