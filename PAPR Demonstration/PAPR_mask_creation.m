load('100MHzLTE.mat');
load('mask.mat');
x = (waveform(1:50000))';
mask_factor = 1;
mask = mask_factor*mask;

mask_db=mag2db(abs(mask));
X=abs(fftshift(fft(x)));
X_db=mag2db(X);

mask_db_opt=mask_db;
%mask_db_opt(2000:4000)=-20;
mask_db_opt(1:2000)=-35;
mask_db_opt(48000:50000)=-35;

mask_opt=(1/20)*exp(mask_db_opt);

figure();
plot(mask_db,'LineWidth',1.5);
hold on;
plot(X_db);
hold on;
plot(mask_db_opt,'LineWidth',1.5);
hold off;
legend('mask','X','mask_opt');

mask_opt=db2mag(mask_db_opt);

%only if mask_opt is good - save it.
%mask=mask_opt;
%save('our_mask_100MHzLTE.mat','mask');