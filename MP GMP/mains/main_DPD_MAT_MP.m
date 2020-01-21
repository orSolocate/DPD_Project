clc; clear; close all; 

%% General User inputs
amp_data          = 'simrfV2_powamp_dpd_data.mat';
signal            = 'testsignal.mat';
model_run_period  = 50000 ;                                %num of samples to run in the model
start_pos_sig     = 1;
end_pos_sig       = start_pos_sig+model_run_period-1;

%% Model User inputs
miu_MP    = 0.005;
mem_depth = 4 ;                                           %M in the MP model
mem_deg   = 7 ;                                           %K in the MP model

%% Loads
load(signal);                                              %load signal
load(amp_data);                                            %load amp data
x               = x(start_pos_sig:end_pos_sig)./(2*norm(x,2));
AMP_coef_Matrix = Get_coef_MP(inDataPA, outDataPA, mem_deg, mem_depth);
G               = mean(abs(outDataPA))/mean(abs(inDataPA));%abs(mean(y./x));
y               = [zeros(mem_depth,1); Get_model_output_MP(AMP_coef_Matrix, x.*G, mem_deg, mem_depth)];

%% Initialize and Model
y    = ifft(fft(y).*exp(-phdiffmeasure(x, y)*1i));
coef = Get_coef_MP(y', x', mem_deg, mem_depth);

%% Itarations pre
k               = 1;
error_decreases = true;
DISPLAY_FACTOR  = 10;
max_iter        = 50;
error_lim       = 1e-5;
iter_error      = ones(1,1).*inf;

%% First iteration
z             = [zeros(mem_depth,1);PD_MP(x, coef, mem_deg, mem_depth)];
z_hat         = [zeros(mem_depth,1);PD_MP(y, coef, mem_deg, mem_depth)];
error         = z - z_hat;
iter_error(k) = norm(error,2);

%% Itarate
while((iter_error(k) > error_lim) && (error_decreases) && (k < max_iter))
    k = k + 1;
    
    delay    = finddelay(error,y);
    coef_hat = Get_coef_MP(y', error', mem_deg, mem_depth);
    coef     = (coef + miu_MP.*coef_hat);  
    
    z               = [zeros(mem_depth,1);PD_MP(x, coef, mem_deg, mem_depth)];
    z_hat           = [zeros(mem_depth,1);PD_MP(y, coef, mem_deg, mem_depth)];
    error           = z - z_hat;
    iter_error(k)   = norm(error,2);
    error_decreases = iter_error(k) <= iter_error(k-1);
    if(rem(k,DISPLAY_FACTOR) == 0)
        disp(['Iteration #: ',num2str(k), ' - Error is: ',num2str(iter_error(k))]);
    end
 
end

%% Plots
figure; plot(iter_error);

z_new   = [zeros(mem_depth,1);PD_MP(x.*G, coef, mem_deg, mem_depth)];

y_MP    = [zeros(mem_depth,1); Get_model_output_MP(AMP_coef_Matrix, z_new.*G, mem_deg, mem_depth)].*G;
y_no_MP = [zeros(mem_depth,1); Get_model_output_MP(AMP_coef_Matrix, x.*G, mem_deg, mem_depth)];

figure; pspectrum(y_MP, 200e6); hold on; pspectrum(y_no_MP, 200e6); legend('y_{MP}','y_{no MP}');
