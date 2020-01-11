clc; clear; close all; 

%% General User inputs
amp_data          = 'simrfV2_powamp_dpd_data.mat';
signal            = 'testsignal.mat';
%signal_verify='100MHzLTE.mat'
model_run_period  = 50000 ;                                %num of samples to run in the model
start_pos_sig     = 1;
end_pos_sig       = start_pos_sig+model_run_period-1;

%% Model User inputs
miu=linspace(-1,1,6);
mem_depth = linspace(0,10,11) ;                                           %M in the MP model
mem_deg   = linspace(1,9,9) ;                                           %K in the MP model

%% Load input
load(signal);                                              %load signal
load(amp_data);                                            %load amp data
x                     = x(start_pos_sig:end_pos_sig)./(2*norm(x,2));
%RMS_in                = 10*log10( norm(x)^2/50/length(x)) + 30 + 20;

%% Optimization
%For this to work you need to be in the directory of this main matlab file
dir=strcat('Optimization/MP_Matlab/',signal(1:end-4));
for u=miu %create a new folder for every MIU (!!!)    
    curr_dir=strcat(dir,'/',num2str(u),'_miu');
    mkdir(curr_dir);
    for i=mem_depth
        for j=mem_deg
            AMP_coef_Matrix  = Get_coef_MP(inDataPA, outDataPA, j, i); 
            %[y, ~, ~, ~] = RFWebLab_PA_meas_v1_1(x, RMS_in);
            y     = [zeros(i,1); Get_model_output_MP(AMP_coef_Matrix, x, j, i)];

            % Initialize and Model
            %WL_delay = finddelay(x, y);
            %if(WL_delay >= 0)
              %  y = y(WL_delay+1:end);
             %   x = x(1:end-WL_delay);
            %elseif(WL_delay < 0)
              %  y = y(1:end-WL_delay);
             %   x = x(WL_delay+1:end);
            %end
            %y_sf = ifft(fft(y).*exp(-phdiffmeasure(x, y)*1i));
            G    = mean(abs(y))/mean(abs(x));%abs(mean(y./x));
            coef = Get_coef_MP(y', x', j, i);

            % Itarations pre
            k          = 1;
            error_decreases  = true;
            DISPLAY_FACTOR=10;
            max_iter   = 50;
            error_lim  = 1e-5;
            iter_error = ones(1,1).*inf;

            % First iteration
            z             = [zeros(i,1);PD_MP(x, coef, j, i)];
            z_hat         = [zeros(i,1);PD_MP(y./G, coef, j, i)];
            error         = z - z_hat;
            iter_error(k) = norm(error,2);

           % Itarate
            while((iter_error(k) > error_lim) && (error_decreases) && (k < max_iter))
                k = k + 1;

                %delay    = finddelay(error,y);
                coef_hat = Get_coef_MP(y', error', j, i);
                coef     = (coef + u.*coef_hat);  

                z             = [zeros(i,1);PD_MP(x, coef, j, i)];
                z_hat         = [zeros(i,1);PD_MP(y./G, coef, j, i)];
                error         = z - z_hat;
                iter_error(k) = norm(error,2);
                error_decreases     = iter_error(k) <= iter_error(k-1);
                if(rem(k,DISPLAY_FACTOR)==0)
                    disp(['Iteration #: ',num2str(k), ' - Error is: ',num2str(iter_error(k))]);
                end

            end

            % Plots
            h(1)=figure;plot(iter_error);

            load(signal);                                              %load signal
            %x=waveform;
            %z = [zeros(i,1);PD_MP(x_new.*G, coef, j, i)];
            z_new = [zeros(i,1);PD_MP(x.*G, coef, j, i)];
            %[y_MP, ~, ~, ~]             = RFWebLab_PA_meas_v1_1(z_new, RMS_in);
             y_MP  = [zeros(i,1); Get_model_output_MP(AMP_coef_Matrix, z_new, j, i)];
             y_MP  = y_MP.*G;
            %[y_no_MP, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_1(x.*G, RMS_in);
             y_no_MP  = [zeros(i,1); Get_model_output_MP(AMP_coef_Matrix, x.*G, j, i)];
            h(2)=figure; pspectrum(y_MP); hold on; pspectrum(y_no_MP); legend('y_{MP}','y_{no MP}');hold off;
            
            save_path=strcat(curr_dir,'/',int2str(i),'_depth_',int2str(j),'_deg_.fig');                       
            savefig(h,save_path,'compact');
            close(h)
        end
    end
end
