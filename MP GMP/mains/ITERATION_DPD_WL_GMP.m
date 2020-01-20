clc; clear; close all; 

%% General User inputs
signal            = 'testsignal.mat';
model_run_period  = 100000 ;                                %num of samples to run in the model
start_pos_sig     = 1;
end_pos_sig       = start_pos_sig+model_run_period-1;

%% Model User inputs
Miu=linspace(0.3,0.6,3);
Ma = linspace(4,5,2) ; Mb = linspace(2,3,2) ; Mc = linspace(2,3,2) ;                                 %memory depth for GMP model
Ka = linspace(8,9,2) ; Kb = linspace(5,6,2) ; Kc = linspace(5,6,2) ;                                 %non-linearity deg for GMP
P = linspace(5,6,2) ;  Q = linspace(5,6,2) ;                                 %cross-terms for GMP

num_iterations=length(Miu)*length(Ma)*length(Mb)*length(Mc)*length(Ka)*length(Kb)*length(Kc)*length(P)*length(Q);
       
%% Loads
load(signal);                                              %load signal
x                     = x(start_pos_sig:end_pos_sig)./(2*norm(x,2));
RMS_in                = 10*log10( norm(x)^2/50/length(x)) + 30 + 20;
[y, ~, ~, ~] = RFWebLab_PA_meas_v1_1(x, RMS_in);

%% Time and Phase Delay fix
WL_delay = finddelay(x, y);
if(WL_delay >= 0)
    y_st = y(WL_delay+1:end);
    x    = x(1:end-WL_delay);
elseif(WL_delay < 0)
    y_st = y(1:end-WL_delay);
    x    = x(WL_delay+1:end);
end
y_sft = ifft(fft(y_st).*exp(-phdiffmeasure(x, y_st)*1i));
G    = mean(abs(y_sft))/mean(abs(x));%abs(mean(y./x));
[y_no_GMP, RMSout, Idc, Vdc] = RFWebLab_PA_meas_v1_1(x.*G, RMS_in);

dir=strcat('..\Optimization\GMP_WL\',signal(1:end-4));

%% Iterations
%For this to work you need to be in the directory of this main matlab file
errors_vector=zeros(num_iterations,1);
idx=0;
for miu=Miu
    curr_dir=strcat(dir,'/',num2str(miu),'_miu');
    mkdir(curr_dir);
    idx=idx+1;
    for ma=Ma
        for mb=Mb
            for mc=Mc
                for ka=Ka
                    for kb=Kb
                        for kc=Kc
                            for p=P
                                for q=Q
                                     orders              = [ma, ka, mb, kb, p, mc, kc, q];
                                     naming=strcat(int2str(ma),'_Ma_',int2str(mb),'_Mb_',int2str(mc),'_Mc_',int2str(ka),'_Ka_',int2str(kb),'_Kb_',int2str(kc),'_Kc_',int2str(p),'_P_',int2str(q),'_Q_');
                                     [iter_error,coef]=PD_iterate(x,y_sft,orders,miu);
                                     z_new = [zeros(first_n_GMP-1,1);PD_GMP(x.*G, coef,orders); zeros(orders(8),1)];
                                     [y_GMP, ~, ~, ~]             = RFWebLab_PA_meas_v1_1(z_new, RMS_in);
                                     h(1)=figure; plot(iter_error);
                                     title("eror for "+naming);
                                     if isempty(y_GMP)==0 && isempty(y_no_GMP)==0
                                        h(2)=figure; 
                                        pspectrum(y_GMP, 200e6); hold on; pspectrum(y_no_GMP, 200e6); legend('y_{MP}','y_{no MP}');
                                        title("results for "+naming);
                                        hold off;
                                     end
                                     save_path=strcat(curr_dir,'/',naming,'.fig');                       
                                     savefig(h,save_path,'compact');
                                     close(h);
                                     %errors_vector(idx)=iter_error;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
disp(["iteration completed",num_iterations]);


function [iter_error,coef]=PD_iterate(x,y_sft,orders,miu)
    coef = Get_coef_GMP(y_sft', x', orders);
    first_n_GMP         = 1+max([orders(1), orders(3)+orders(5), orders(6)-1]);

    k                = 1;
    error_decreases  = true;
    DISPLAY_FACTOR   = 10;
    max_iter         = 50;
    error_lim        = 1e-5;
    iter_error       = ones(1,1).*inf;

    z             = [zeros(first_n_GMP-1,1);PD_GMP(x, coef,orders); zeros(orders(8),1)];
    z_hat         = [zeros(first_n_GMP-1,1);PD_GMP(y_sft, coef,orders); zeros(orders(8),1)];
    error         = z - z_hat;
    iter_error(k) = norm(error,2);

    while((iter_error(k) > error_lim) && (error_decreases) && (k < max_iter))
        k = k + 1;

        y_sft = ifft(fft(y_sft).*exp(-phdiffmeasure(error, y_sft)*1i));
        delay      = finddelay(error,y_sft);
        if(delay >= 0)
            y_sft = y_sft(delay+1:end);
            error = error(1:end-delay);
        elseif(delay < 0)
            y_sft = y_sft(1:end-delay);
            error = error(delay+1:end);
        end

        coef_hat = Get_coef_GMP(y_sft', error', orders);
        coef     = (coef + miu.*coef_hat);  
        y_sft    = [y_sft; zeros(length(x)-length(y_sft),1)];

        z               = [zeros(first_n_GMP-1,1);PD_GMP(x, coef,orders); zeros(orders(8),1)];
        z_hat           = [zeros(first_n_GMP-1,1);PD_GMP(y_sft, coef,orders); zeros(orders(8),1)];
        error           = z - z_hat;
        iter_error(k)   = norm(error,2);
        error_decreases = iter_error(k) <= iter_error(k-1);

        if(rem(k,DISPLAY_FACTOR)==0)
            disp(['Iteration #: ',num2str(k), ' - Error is: ',num2str(iter_error(k))]);
        end
    end
end