clc; clear; close all; 

%% General User inputs
signal            = 'testsignal.mat';
model_run_period  = 100000 ;                                %num of samples to run in the model
start_pos_sig     = 1;
end_pos_sig       = start_pos_sig+model_run_period-1;

%% Model User inputs
miu_GMP    = -1;
Ma = 8 ; Mb = 8 ; Mc = 8 ;                                 %memory depth for GMP model
Ka = 9 ; Kb = 5 ; Kc = 5 ;                                 %non-linearity deg for GMP
          P = 3 ;  Q = 3 ;                                 %cross-terms for GMP
          
%% Loads
load(signal);                                              %load signal
x                     = x(start_pos_sig:end_pos_sig)./(2*norm(x,2));
RMS_in                = 10*log10( norm(x)^2/50/length(x)) + 30 + 20;
[y, ~, ~, ~] = RFWebLab_PA_meas_v1_1(x, RMS_in);

%% Initialize and Model
orders              = [Ma, Ka, Mb, Kb, P, Mc, Kc, Q];
first_n_GMP         = 1+max([orders(1), orders(3)+orders(5), orders(6)-1]);
WL_delay = finddelay(x, y);
if(WL_delay >= 0)
    y = y(WL_delay+1:end);
    x = x(1:end-WL_delay);
elseif(WL_delay < 0)
    y = y(1:end-WL_delay);
    x = x(WL_delay+1:end);
end
y_sf = ifft(fft(y).*exp(-phdiffmeasure(x, y)*1i));
G    = mean(abs(y_sf))/mean(abs(x));%abs(mean(y./x));
coef = Get_coef_GMP(y_sf', x', orders);

%% Itarations pre
k                = 1;
error_decreases  = true;
DISPLAY_FACTOR   = 10;
max_iter         = 50;
error_lim        = 1e-5;
iter_error       = ones(1,1).*inf;

%% First iteration
z             = [zeros(first_n_GMP-1,1);PD_GMP(x, coef,orders); zeros(orders(8),1)];
z_hat         = [zeros(first_n_GMP-1,1);PD_GMP(y_sf, coef,orders); zeros(orders(8),1)];
error         = z - z_hat;
iter_error(k) = norm(error,2);

%% Itarate
while((iter_error(k) > error_lim) && (error_decreases) && (k < max_iter))
    k = k + 1;
    
    delay      = finddelay(error,y_sf);
    if(WL_delay >= 0)
        y_sf  = y_sf(WL_delay+1:end);
        error = error(1:end-WL_delay);
    elseif(WL_delay < 0)
        y_sf  = y_sf(1:end-WL_delay);
        error = error(WL_delay+1:end);
    end
    
    coef_hat = Get_coef_GMP(y_sf', error', orders);
    coef     = (coef + miu_GMP.*coef_hat);  
    y_sf     = [y_sf; zeros(length(x)-length(y_sf),1)];
    
    z               = [zeros(first_n_GMP-1,1);PD_GMP(x, coef,orders); zeros(orders(8),1)];
    z_hat           = [zeros(first_n_GMP-1,1);PD_GMP(y_sf, coef,orders); zeros(orders(8),1)];
    error           = z - z_hat;
    iter_error(k)   = norm(error,2);
    error_decreases = iter_error(k) <= iter_error(k-1);
    
    if(rem(k,DISPLAY_FACTOR)==0)
        disp(['Iteration #: ',num2str(k), ' - Error is: ',num2str(iter_error(k))]);
    end
 
end

%% Plots
figure; plot(iter_error);

z_new = [zeros(first_n_GMP-1,1);PD_GMP(x.*G, coef,orders); zeros(orders(8),1)];

[y_GMP, ~, ~, ~]             = RFWebLab_PA_meas_v1_1(z_new, RMS_in);
[y_no_GMP, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_1(x, RMS_in);
 
figure; pspectrum(y_GMP); hold on; pspectrum(y_no_GMP); legend('y_{GMP}','y_{no GMP}');
