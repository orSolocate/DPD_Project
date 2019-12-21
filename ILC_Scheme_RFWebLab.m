function [x_opt, error_vec_plot] = ILC_Scheme_RFWebLab(x,y,orders,N)

% Error requirement

err_req = 1e-9;
error_vec_plot = [];
sigma=0.001;
% Block to get PA model

WL_delay               = finddelay(x, y);

if(WL_delay >=0)
    avg_gain_weblab           = abs(mean(y(WL_delay+1:end)./x(1:end-WL_delay)));
    PD_coef_Matrix_WL  = Get_coef_MP(y(WL_delay+1:end)', x(1:end-WL_delay)', orders(1), orders(2));
end
if(WL_delay < 0)
    avg_gain_weblab           = abs(mean(y(1:end-WL_delay)./x(WL_delay+1:end)));
    PD_coef_Matrix_WL  = Get_coef_MP(y(1:end-WL_delay)', x(WL_delay+1:end)',orders(1), orders(2));
end
y_d_WL                 = x.*avg_gain_weblab;
x_opt                  = zeros(length(y_d_WL),1);

%numDataPts = length(y_d_WL);
%halfDataPts = round(numDataPts/2); % Half is relevant for simulation purposes only!

%coefMat = Get_Coefficients_Matrix(inDataPA(1:halfDataPts), outDataPA(1:halfDataPts), memLen, degLen);

% error calculation
n = orders(2)+1;
while (n<length(y_d_WL))                                         %run the system as live
    x_opt(n) = PD_MP(y_d_WL(n-orders(2):n), PD_coef_Matrix_WL, orders(1),orders(2));
    n                  = n+1;
end

u=y_d_WL./avg_gain_weblab;

y_d_WL=y_d_WL(WL_delay+1:end) %get rid of delay

y =-1* Get_model_output_MP(PD_coef_Matrix_WL, u, orders(1), orders(2))
%noise_i = normrnd(zeros(1,length(y)), sigma*ones(1,length(y)));
%noise_q = normrnd(zeros(1,length(y)), sigma*ones(1,length(y)));
%y = y + noise_i + 1i * noise_q;
error_vec = y_d_WL(orders(1):length(y_d_WL)) - y(orders(1):length(y_d_WL));
%error_vec=y_d_WL(WL_delay+1:end)-y
error = norm(error_vec, 2);
error_vec_plot = [error_vec_plot error];

gamma = 2/(2*avg_gain_weblab);

for ii=1:N
    if (error <= err_req)
        break;
    end 
    u(orders(1):length(y_d_WL)) = u(orders(1):length(y_d_WL)) + gamma * error_vec;
    % error calculation
    max(abs(u))
    y = -1*Get_model_output_MP(PD_coef_Matrix_WL, u, orders(1), orders(2))
    %u=u./avg_gain_weblab;
    %[y, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(u);
    pspectrum(y)
    error_vec = y_d_WL(orders(1):length(y_d_WL)) -y(orders(1):length(y_d_WL));
    %error_vect=y_d_WL(WL_delay+1:end)-y
    error = norm(error_vec, 2);
    error_vec_plot = [error_vec_plot error];
end

end