%% initializations

load('100MHzLTE.mat');
x_hat = waveform(1:10000);
x_hat_norm = norm(x_hat,2);
[y_hat, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(x_hat./(10*x_hat_norm));
WL_delay               = finddelay(x_hat, y_hat);

if(WL_delay >=0)
    g_avg           = abs(mean(y_hat(WL_delay+1:end)./x_hat(1:end-WL_delay)));
end
if(WL_delay < 0)
    g_avg           = abs(mean(y_hat(1:end-WL_delay)./x_hat(WL_delay+1:end)));
end
yd = (waveform(1:10000)) * g_avg;

buffer = [];
predistorter_iter = 3000;
u_star_vec = zeros(1, predistorter_iter);

%% ILC_Scheme

[u, error_vec_plot] = ILC_Scheme(yd, g_avg, x_hat, y_hat, x_hat_norm );
figure; 
plot(1:length(error_vec_plot), error_vec_plot);
xlabel('Iterations');
ylabel('Error');
title('Error vs. Iteration Number');
grid;
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','GridAlpha'),'GridAlpha',1);

%% System Identification

coefMat = System_Identification(yd, u);

%% Predistorter

yd = (waveform(55000:65000))' * g_avg; %Changing yd to apply general case (online input)

for ii = 1:predistorter_iter
    [u_star, buffer] = Predistorter(yd(ii)/g_avg, buffer, coefMat, g_avg);
    u_star_vec(ii) = u_star;
end


%% Comparison

numDataPts = length(x_hat);
memLen = 3;
degLen = 7;

halfDataPts = round(numDataPts/2);

amplifier = Get_Coefficients_Matrix(x_hat(1:halfDataPts), y_hat(1:halfDataPts), memLen, degLen);


[y1, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(yd(1:length(u_star_vec))/g_avg);

[y2, RMSout, Idc, Vdc]  = RFWebLab_PA_meas_v1_1(u_star_vec);

figure;
plot(abs(yd(1:length(u_star_vec))/g_avg), abs(y1), 'o');
hold all;
plot(abs(yd(1:length(u_star_vec))/g_avg), abs(y2), 'o');
legend('Without ILC-DPD','With ILC-DPD')
xlabel('abs(input)');
ylabel('abs(output)');
title('Linearity Comparison');
grid;
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','GridAlpha'),'GridAlpha',1);
