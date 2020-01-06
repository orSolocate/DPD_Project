function x_CFR = PAPR(x, Th, mask, Amax, Smax) 
% Amax should be a parameter, constant or we can try few values and compare
% the results
%EVM Threshold (Th) and Spectral mask (P[k]) also should be constants or function parameters
%I didnt use fftshift, not sure if actually needed

L = 4; %According to the paper, L is typically >=4
N = length(x);

% Adjust Spectral Mask to Xk_L
P = zeros(1,L*N);
P(1:N) = fftshift(mask);
P(N+1:L*N) = mask(1);

% The upsampled x is level is lower by L (due to fft/ifft formula)
Amax = Amax/L;

%defining series
Xk = fft(x);
Xk_CFR = zeros(1,L*N);

%Upsampling
Xk_L = zeros(1,L*N);
Xk_L(1:N) = Xk;
xn_L = L*ifft(Xk_L);

%clipping
xn_L_clipped = xn_L;
xn_L_clipped(abs(xn_L)>Amax) = Amax;
Xk_L_clipped = fft(xn_L_clipped);

%In-band processing
Err_k = Xk_L_clipped(1:N) - Xk;
if EVM(Err_k, Smax)<Th
    Xk_CFR(1:N) = Xk_L_clipped(1:N);
else
    [sorted_Err_k, Sorted_Index] = sort(abs(Err_k));
    m=1;
    while(rms(sorted_Err_k(1:m)) < Th*Smax)
        m = m+1;
    end
    M = m-1;
    M_Set = Sorted_Index(1:M);
    I_sub_M = Sorted_Index(M+1:N);
    Xk_CFR(M_Set) = Xk_L_clipped(M_Set);
    Xk_CFR(I_sub_M) = Xk(I_sub_M) + Th*Smax*exp(1j*angle(Err_k(I_sub_M)));  
end

%Out-Of-Band processing
for k = N+1:L*N
    if (abs(Xk_L_clipped(k))^2 <= abs(P(k))^2)
       Xk_CFR(k) = Xk_L_clipped(k);
    else
        Xk_CFR(k) = P(k)*exp(1j*angle(Xk_L_clipped(k)));
    end
end

x_CFR = downsample(ifft(Xk_CFR),L); %% when we upsmapled the fft terms decreased by L

end