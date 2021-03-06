clc; clear; close all; 

%% General User inputs
signal            = 'testsignal.mat';
model_run_period  = 100000 ;                                %num of samples to run in the model
start_pos_sig     = 1;
end_pos_sig       = start_pos_sig+model_run_period-1;

%% Model User inputs
miu_MP    = 0.01;
mem_depth = 5 ;                                           %M in the MP model
mem_deg   = 7 ;                                           %K in the MP model

noise_var = 0.001;

%% Loads
load(signal);                                              %load signal
x            = x(start_pos_sig:end_pos_sig)./(2*norm(x,2));

%% PAPR stuff
% load('mask.mat');
% Amax_factor = 1; %5
% mask_factor = 1; %5
% Smax = sqrt(98/42);
% Th = 0.001 * EVM(x, Smax); %0.08
% mask = mask_factor*mask;
% Amax = Amax_factor*mean(abs(x));
% x = (PAPR(x', Th, mask, Amax, Smax))'; % works for signals length of 50k

%% Initialize and Model
RMS_in       = 10*log10( norm(x)^2/50/length(x)) + 30 + 20;
[y, ~, ~, ~] = RFWebLab_PA_meas_v1_1(x, RMS_in);

WL_delay = finddelay(x, y);
if(WL_delay >= 0)
    y_st = y(WL_delay+1:end);
    x    = x(1:end-WL_delay);
elseif(WL_delay < 0)
    y_st = y(1:end+WL_delay);
    x    = x(1-WL_delay:end);
end

y_sft = ifft(fft(y_st).*exp(-phdiffmeasure(x, y_st)*1i));
G     = mean(abs(y_sft))/mean(abs(x));%abs(mean(y./x));
noise = wgn(length(y_sft),1,10*log10(noise_var),'complex');
coef  = Get_coef_MP((y_sft + noise)', x', mem_deg, mem_depth);

%% Itarations pre
k                = 1;
error_decreases  = true;
DISPLAY_FACTOR   = 10;
max_iter         = 50;
error_lim        = 1e-5;
iter_error       = ones(1,1).*inf;

%% First iteration
z             = [zeros(mem_depth,1);PD_MP(x, coef, mem_deg, mem_depth)];
z_hat         = [zeros(mem_depth,1);PD_MP(y_sft, coef, mem_deg, mem_depth)];
error         = z - z_hat;
iter_error(k) = norm(error,2);

%% Itarate
while((iter_error(k) > error_lim) && (error_decreases) && (k < max_iter))
    k = k + 1;
    
    y_sft = ifft(fft(y_sft).*exp(-phdiffmeasure(error, y_sft)*1i));
    delay = finddelay(error,y_sft);
    if(delay >= 0)
        coef_hat = Get_coef_MP((y_sft(1+delay:end))', error(1:end-delay)', mem_deg, mem_depth);
    elseif(delay < 0)
        coef_hat = Get_coef_MP((y_sft(1:end+delay))', error(1-delay:end)', mem_deg, mem_depth);
    end 
    coef = (coef + miu_MP.*coef_hat);  
    
    z               = [zeros(mem_depth,1);PD_MP(x, coef, mem_deg, mem_depth)];
    z_hat           = [zeros(mem_depth,1);PD_MP(y_sft, coef, mem_deg, mem_depth)];
    error           = z - z_hat;
    iter_error(k)   = norm(error,2);
    error_decreases = iter_error(k) <= iter_error(k-1);
    
    if(rem(k,DISPLAY_FACTOR) == 0)
        disp(['Iteration #',num2str(k), ' - Error is: ',num2str(iter_error(k))]);
    end
 
end

%% Plots
figure; plot(iter_error);

z_new = [zeros(mem_depth,1);PD_MP(x.*G, coef, mem_deg, mem_depth)]; %x.*G*0.9

[y_MP, ~, ~, ~]             = RFWebLab_PA_meas_v1_1(z_new, RMS_in);
[y_no_MP, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_1(x, RMS_in);
 
figure; pspectrum(y_MP, 200e6); hold on; pspectrum(y_no_MP, 200e6); legend('y_{MP}','y_{no MP}');
