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

%% Loads
load(signal);                                              %load signal
x            = x(start_pos_sig:end_pos_sig)./(2*norm(x,2));
RMS_in       = 10*log10( norm(x)^2/50/length(x)) + 30 + 20;
[y, ~, ~, ~] = RFWebLab_PA_meas_v1_1(x, RMS_in);

%% Initialize and Model
WL_delay = finddelay(x, y);
if(WL_delay >= 0)
    y_st = y(WL_delay+1:end);
    x    = x(1:end-WL_delay);
elseif(WL_delay < 0)
    y_st = y(1:end+WL_delay);
    x    = x(1-WL_delay:end);
end

y_sft = ifft(fft(y_st).*exp(-phdiffmeasure(x, y_st)*1i));

%% MATLAB modeling

estimator = comm.DPDCoefficientEstimator(...
    'DesiredAmplitudeGaindB',10, 'PolynomialType','Memory polynomial','Degree',mem_deg,'MemoryDepth',mem_depth,'Algorithm','Least squares');
coef = estimator(x,y_sft);
dpd = comm.DPD('PolynomialType','Memory polynomial','Coefficients',coef);
z = dpd(x.*300);

[y_MP, ~, ~, ~]             = RFWebLab_PA_meas_v1_1(z, RMS_in);
figure; pspectrum(y_MP, 200e6); hold on; pspectrum(y, 200e6); legend('y_{MP}','y_{no MP}');
