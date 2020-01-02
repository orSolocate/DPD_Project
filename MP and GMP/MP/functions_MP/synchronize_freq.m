function [ x2_sync, phase_coeff, phase_offset ] = synchronize_freq( x1, x2, mode, varargin)
% %
%
%   Synchronize two complex datastreams by adjusting correcting the linear
%   phase error
%   Input parameters:   x1, x2 : complex vectors of sam size containing the
%                               data
%                       mode:   'v' - verbose mode
%                               'q' - quiet mode

% convert both input streams to freq domain. flip both halfs of the output
% vector to put dc frequency to the center

%
% linear phase coeff: -1i*2*pi/n_samples<phase_coeff<1i*2*pi/n_samples
%                       (signal already synchronized in time domain at full
%                       samples)
%
%% Setup
POWER_THRESHOLD = -14;
DC_CLEAR        = 1000;%round(numel(x1)*(1e-3)); %0.1% around dc carrier

% for var_arg_sel=1:2:nargin-3
%     if strcmp(varargin{var_arg_sel},            'POWER_THRESHOLD')
%         POWER_THRESHOLD = varargin{var_arg_sel+1};
%     elseif strcmp(varargin{var_arg_sel},        'DC_CLEAR')
%         DC_CLEAR = varargin{var_arg_sel+1};
%     end
% end

%% Default output
x2_sync         = 0;
phase_coeff     = 0;
%% compatibility issues
phase_offset    = 0;
%

if size(x1,1)~=size(x2,1) || numel(x1)~=numel(x2)
    disp('Error: Please make sure x1 and x2 are of the same dimension.')   
    return
end

TRANSPOSE = false;
if size(x1,1)>size(x1,2)
    TRANSPOSE = true;
    x1 = x1.';
    x2 = x2.';
end

%%
X1              = fftshift(fft(x1));        % Transform reference signal
X1              = X1/max(abs(X1));          % Normalize reference signal

X2              = fftshift(fft(x2));        % Transform measured signal
scaling_factor  = max(abs(X2));
X2              = X2/scaling_factor;        % Normalize measured signal
NFFT            = numel(X1);

select_vect     = 20*log10(abs(X1)) > POWER_THRESHOLD;      % select only high power frequencies (-> used channels)
fu_neg          = round(numel(select_vect)/2)-DC_CLEAR;
fu_pos          = round(numel(select_vect)/2)+DC_CLEAR;
if fu_neg>0 && fu_pos<numel(select_vect)
    select_vect(fu_neg:fu_pos) = false;                     % clear area around dc. there are many errors due to carrier leakage
else
    disp('Warning: x1 and x2 are too short to apply DC-clearance.')
    return
end

if sum(select_vect)<2
    disp('Error: no usable signal levels detected. Synchronization deactivated.')
    x2_sync = x2;
    if TRANSPOSE
        x2_sync = x2_sync.';
    end
    return
end

%% UPDATE: 2016-11-28
%% computate lower and upper sideband seperated

LOW         = 1;
HIGH        = 2;
phase_coeff  = zeros(size(X1));
phase_error = zeros(size(X1));

sb_range = 1:numel(X1);%2017_11_03
%for sb_sel=LOW:HIGH
%    if sb_sel == LOW
%        sb_range = 1:numel(X1)/2-1;
%    else
%        sb_range = numel(X1)/2+1:numel(X1);
%    end
    select_vect_tmp             = logical(zeros(size(select_vect)));
    select_vect_tmp(sb_range)   = select_vect(sb_range);
    select_index                = 1:numel(select_vect_tmp);
    
    % compensate linear phase error
    phase_error(sb_range)           = angle(X1(sb_range)./X2(sb_range));
    phase_error(select_vect_tmp)    = unwrap(phase_error(select_vect_tmp));
    
    P                       = polyfit(select_index(select_vect_tmp),phase_error(select_vect_tmp),1);
    phase_coeff_tmp          = polyval(P,sb_range);
    X2(sb_range)            = X2(sb_range).*exp(1i*phase_coeff_tmp);
    
    phase_coeff(sb_range)	= phase_coeff_tmp;
%end
x2_sync = ifft(fftshift(X2*scaling_factor));
if TRANSPOSE
    x2_sync     = x2_sync.';
    phase_coeff = phase_coeff.';

end
    
if strcmp(mode,'v')
    figure(1001)
    % plot the differences in the absolute values to show the amplitude error. x2 will be scaled.
    subplot(4,1,1)
    freq = -0.5:1/(NFFT-1):0.5;
    plot(freq,20*log10(abs(X1)),freq,20*log10(abs(X2)),freq,50*(select_vect-1))
    axis([-0.5 0.5 -60 10])
    ylabel('PSD (dBx/Hz)')
    legend('|original|','|measured| (scaled)','mask')
    grid on
    
    % plot the unwrapped phase difference between x1 and x2
    subplot(4,1,2)
    plot(freq,phase_error)
    legend('phase error (full band)')
    ylabel('Phase error (rad)')
    grid on
    
    % plot the unwrapped phase difference between x1 and x2
    subplot(4,1,3)
    plot(freq,phase_coeff)   % Zero at 0*fs
    legend('phase\_error compensation (full band)')
    ylabel('Phase error (rad)')
    grid on
    
    % plot the phase difference between x1 and x2 after correction of
    % linear error
    subplot(4,1,4)
    plot(freq,angle(X1./X2))
    legend('phase error after correction (full band)')
    ylabel('Phase error (rad)')
    grid on
    %legend('arg\{X1/X2\_sync\}','RMS')
end