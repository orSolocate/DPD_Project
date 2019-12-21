function [x_opt, error_vec_plot] = ILC_Scheme_RFWebLab(x,y,orders,N)

% Error requirement

err_req = 1e-9;
error_vec_plot = [];

% Block to get PA model


weblab_delay           = finddelay(x, y);
if(weblab_delay >=0)
    avg_gain_weblab    = abs(mean(y(weblab_delay+1:end)./x(1:end-weblab_delay)));
    PD_coef_Matrix_WL  = Get_coef_MP(y(weblab_delay+1:end)', x(1:end-weblab_delay)', orders(1), orders(2));
end
if(weblab_delay < 0)
    avg_gain_weblab    = abs(mean(y(1:end-weblab_delay)./x(weblab_delay+1:end)));
    PD_coef_Matrix_WL  = Get_coef_MP(y(1:end-weblab_delay)', x(weblab_delay+1:end)', orders(1), orders(2));
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

u=x./avg_gain_weblab;

[y, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(u);

error_vec = y_d_WL(orders(1):length(y_d_WL)) - y(orders(1):length(y_d_WL));
error = norm(error_vec, 2);
error_vec_plot = [error_vec_plot error];

gamma = 2/(2*avg_gain_weblab);

for ii=1:N
    if (error <= err_req)
        break;
    end 
    u(orders(1):length(y_d_WL)) = u(orders(1):length(y_d_WL)) + gamma * error_vec;
    % error calculation
    %max(abs(u))
    u=u./avg_gain_weblab;
    [y, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(u);
    pspectrum(y)
    error_vec = y_d_WL(orders(1):length(y_d_WL)) - y(orders(1):length(y_d_WL));
    error = norm(error_vec, 2);
    error_vec_plot = [error_vec_plot error];
end

end