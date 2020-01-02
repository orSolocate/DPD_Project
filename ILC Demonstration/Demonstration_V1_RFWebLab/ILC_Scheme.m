function [u, error_vec_plot] = ILC_Scheme(yd, g_avg, x_hat, y_hat, avg_power  )

% Error requirement
iterations_num =10;
err_req = 1e-9;
error_vec_plot = [];
sigma = 0.001;
ii = 0;

% Block to get PA model

numDataPts = length(x_hat);
memLen = 3;
degLen = 7;
halfDataPts = round(numDataPts/2); % Half is relevant for simulation purposes only!
%x_hat=x_hat./x_hat_norm;
coefMat = Get_Coefficients_Matrix(x_hat(1:halfDataPts), y_hat(1:halfDataPts), memLen, degLen);

% "k=1"
u = (yd./g_avg);%./x_hat_norm;

% error calculation
%[y, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(u);%./(x_hat_norm));
[y, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(u,avg_power);
%y = -1*Get_PA_Output(u, coefMat)';
noise_i = normrnd(zeros(1,length(y)), sigma*ones(1,length(y)));
noise_q = normrnd(zeros(1,length(y)), sigma*ones(1,length(y)));
y = y + (noise_i)' + 1i * (noise_q');
error_vec = yd(memLen:length(yd)) - y(memLen:length(yd));
error = norm(error_vec, 2);
error_vec_plot = [error_vec_plot error];

gamma = 2/(2*g_avg);

while (error > err_req)
    ii = ii + 1;
    if (ii == iterations_num)
        break;
    end
    % "iteration k"
    u(memLen:length(yd)) = u(memLen:length(yd)) + gamma * error_vec;
    % error calculation
    max(abs(u))
    %y = -1*Get_PA_Output(u, coefMat)';
    %[y, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(u./(x_hat_norm));
    [y, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(u,avg_power);

    pspectrum(y)
    noise_i = normrnd(zeros(1,length(y)), sigma*ones(1,length(y)));
    noise_q = normrnd(zeros(1,length(y)), sigma*ones(1,length(y)));
    y = y;%% + noise_i + j * noise_q;
    error_vec = yd(memLen:length(yd)) - y(memLen:length(yd));
    error = norm(error_vec, 2);
    error_vec_plot = [error_vec_plot error];
end

end