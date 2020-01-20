clc; clear; close all; 

%% General User inputs
signal            = 'testsignal.mat';
model_run_period  = 100000 ;                                %num of samples to run in the model
start_pos_sig     = 1;
end_pos_sig       = start_pos_sig+model_run_period-1;

%% Model User inputs
Miu=linspace(0.45,0.55,2);
Mem_depth = linspace(5,6,2) ;                                           %M in the MP model
Mem_deg   = linspace(3,4,2) ;                                            %K in the MP model

num_iterations=length(Mem_deg)*length(Mem_depth)*length(Miu);

%% Loads
load(signal);                                              %load signal
x                     = x(start_pos_sig:end_pos_sig)./(2*norm(x,2));
RMS_in                = 10*log10( norm(x)^2/50/length(x)) + 30 + 20;
[y, ~, ~, ~] = RFWebLab_PA_meas_v1_1(x, RMS_in);

% PAPR stuff

%% Time and Phase Delay fix
WL_delay = finddelay(x, y);
if(WL_delay >= 0)
    y_st = y(WL_delay+1:end);
    x = x(1:end-WL_delay);
elseif(WL_delay < 0)
    y_st = y(1:end+WL_delay);
    x = x(1-WL_delay:end);
end

y_sft = ifft(fft(y_st).*exp(-phdiffmeasure(x, y_st)*1i));
G    = mean(abs(y_sft))/mean(abs(x));%abs(mean(y./x));
[y_no_MP, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_1(x.*G, RMS_in);

dir=strcat('..\Optimization/MP_WL/',signal(1:end-4));

%% Iterations
%For this to work you need to be in the directory of this main matlab file

idx=0;
for miu=Miu
    curr_dir=strcat(dir,'/',num2str(miu),'_miu');
    mkdir(curr_dir);
    idx=idx+1;
    for mem_depth=Mem_depth
        for mem_deg=Mem_deg
            naming=strcat(num2str(mem_depth)," depth ",num2str(mem_deg)," deg");
            
            [iter_error,coef]=PD_iterate(x,y_sft,mem_depth,mem_deg,miu);
            
            z_new = [zeros(mem_depth,1);PD_MP(x.*G, coef, mem_deg, mem_depth)];
            [y_MP, ~, ~, ~] = RFWebLab_PA_meas_v1_1(z_new, RMS_in);
             
            h(1)=figure;plot(iter_error);
            title("eror for "+naming);
            if isempty(y_MP)==0 && isempty(y_no_MP)==0
                h(2)=figure; 
                pspectrum(y_MP, 200e6); hold on; pspectrum(y_no_MP, 200e6); legend('y_{MP}','y_{no MP}');
                title("results for "+naming);
                hold off;
            end
            save_path=strcat(curr_dir,'/',naming,'.fig');                       
            savefig(h,save_path,'compact');
            close(h)
        end
    end
end
disp(["iteration completed",num_iterations]);


function [iter_error,coef]=PD_iterate(x,y_sft,mem_depth,mem_deg,miu)
    coef = Get_coef_MP(y_sft', x', mem_deg, mem_depth);
    k                = 1;
    error_decreases  = true;
    DISPLAY_FACTOR   = 10;
    max_iter         = 50;
    error_lim        = 1e-5;
    iter_error       = ones(1,1).*inf;

    % First iteration
    z             = [zeros(mem_depth,1);PD_MP(x, coef, mem_deg, mem_depth)];
    z_hat         = [zeros(mem_depth,1);PD_MP(y_sft, coef, mem_deg, mem_depth)];
    error         = z - z_hat;
    iter_error(k) = norm(error,2);

    % Itarate
    while((iter_error(k) > error_lim) && (error_decreases) && (k < max_iter))
        k = k + 1;

        delay    = finddelay(error,y_sft);
        if(delay >= 0)
            coef_hat = Get_coef_MP((y_sft(1+delay:end))', error(1:end-delay)', mem_deg, mem_depth);
        elseif(delay < 0)
            coef_hat = Get_coef_MP((y_sft(1:end+delay))', error(1-delay:end)', mem_deg, mem_depth);
        end 
        coef     = (coef + miu.*coef_hat);  

        z               = [zeros(mem_depth,1);PD_MP(x, coef, mem_deg, mem_depth)];
        z_hat           = [zeros(mem_depth,1);PD_MP(y_sft, coef, mem_deg, mem_depth)];
        error           = z - z_hat;
        iter_error(k)   = norm(error,2);
        error_decreases = iter_error(k) <= iter_error(k-1);

        if(rem(k,DISPLAY_FACTOR) == 0)
            disp(['Iteration #',num2str(k), ' - Error is: ',num2str(iter_error(k))]);
        end
    end
end