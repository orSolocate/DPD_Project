% Configure the set of NDLRB values to describe the carriers to be
% aggregated
NDLRB = [100 100 100 100 100];

% Establish the number of component carriers
numCC = numel(NDLRB);

if numCC<2
    error('Please specify more than one CC bandwidth value.');
end

%%
% Configure the number of subframes to generate
numSubframes = 10;

% CC configuration
enb = cell(1,numCC);
for i = 1:numCC
    switch NDLRB(i)
        case 6
            enb{i} = lteRMCDL('R.4');
        case 15
            enb{i} = lteRMCDL('R.5');
        case 25
            enb{i} = lteRMCDL('R.6');
        case 50
            enb{i} = lteRMCDL('R.7');
        case 75
            enb{i} = lteRMCDL('R.8');
        case 100
            enb{i} = lteRMCDL('R.9');
        otherwise
            fprintf('Not a valid number of resource blocks: %d\n',...
                NDLRB(i));
            return
    end
    enb{i}.NDLRB = NDLRB(i);
    enb{i}.Bandwidth = hNRBToBandwidth(NDLRB(i));
    enb{i}.TotSubframes = numSubframes;
    enb{i}.PDSCH.PRBSet = (0:enb{i}.NDLRB-1).';
    enb{i}.PDSCH.RVSeq = 0;
    enb{i}.NCellID = 10;
end

%%
cec = struct;                        % Channel estimation config structure
cec.PilotAverage = 'UserDefined';    % Type of pilot symbol averaging
cec.FreqWindow = 15;                 % Frequency window size
cec.TimeWindow = 15;                 % Time window size
cec.InterpType = 'Cubic';            % 2D interpolation type
cec.InterpWindow = 'Centered';       % Interpolation window type
cec.InterpWinSize = 1;               % Interpolation window size
%%
% Define Delta_f1 parameter in MHz for nominal guard band and frequency
% offset calculations. In the downlink Delta_f1 is the subcarrier spacing
% (TS 36.101, Table 5.6A-1)
Delta_f1 = 0.015; % in MHz, LTE subcarrier spacing in MHz
maxBW = hNRBToBandwidth(max(NDLRB));


F_c = zeros(1,numCC);

ccSpacing = zeros(1,numCC-1); % CC spacing

%  Calculate nominal guard band, TS 36.101 5.6A-1
%Or: formula for aggregated transmission bandwith 25 to 100:
nominalGuardBand = 0.05*maxBW-0.5*Delta_f1;

% Initially assume lower carrier frequency is at baseband
F_c(1) = 0;

% Calculate CC spacing and carrier values
for k = 2:numCC
    ccSpacing(k-1) = hCarrierAggregationChannelSpacing( ...
        enb{k-1}.Bandwidth, enb{k}.Bandwidth);
    F_c(k) = F_c(k-1) + ccSpacing(k-1);
end

% Calculate lower and higher frequency offsets, TS 36.101 5.6A
F_offset_low = (0.18*NDLRB(1)+Delta_f1)/2 + nominalGuardBand;
F_offset_high = (0.18*NDLRB(end)+Delta_f1)/2 + nominalGuardBand;

% Calculate lower and higher frequency edges, TS 36.101 5.6A
F_edge_low = F_c(1) - F_offset_low;
F_edge_high = F_c(end) + F_offset_high;

% Calculate aggregated channel bandwidth, TS 36.101 5.6A
BW_channel_CA = F_edge_high - F_edge_low;
fprintf('BW_channel_CA: %0.4f MHz\n',BW_channel_CA);

% Calculate shift to center baseband
shiftToCenter = -1*(BW_channel_CA/2 + F_edge_low);

% Aggregated bandwidth centered at baseband
F_c = F_c + shiftToCenter;
F_edge_low = F_c(1) - F_offset_low;
F_edge_high = F_c(end) + F_offset_high;

% Display frequency band edges
fprintf('F_edge_low:  %0.4f MHz\n',F_edge_low);
fprintf('F_edge_high: %0.4f MHz\n',F_edge_high);
fprintf('F_offset_low:  %0.4f MHz\n',F_offset_low);
fprintf('F_offset_high: %0.4f MHz\n',F_offset_high);

% Display carrier frequencies
fprintf('\n');
for i = 1:numCC
    fprintf('Component Carrier %d:\n',i);
    fprintf('   Fc: %0.4f MHz\n', F_c(i));
end

%%
% Bandwidth utilization of 85%
bwfraction = 0.85;

% Calculate sampling rates of the component carriers
CCSR = zeros(1,numCC);
for i = 1:numCC
    info = lteOFDMInfo(enb{i});
    CCSR(i) = info.SamplingRate;
end

% Calculate overall sampling rate for the aggregated signal

% Calculate the oversampling ratio (for the largest BW CC) to make sure the
% signal occupies 85% (bwfraction) of the total bandwidth
OSR = (BW_channel_CA/bwfraction)/(max(CCSR)/1e6);
% To simplify the resampling operation choose an oversampling ratio which
% is a power of 2: calculate the next power of two above OSR
OSR = 2^ceil(log2(OSR));
SR = OSR*max(CCSR);
fprintf('\nOutput sample rate: %0.4f Ms/s\n\n',SR/1e6);

% Calculate individual oversampling factors for the component carriers
OSRs = SR./CCSR;

%%
% Generate component carriers
tx = cell(1,numCC);
for i = 1:numCC
    tx{i} = lteRMCDLTool(enb{i},randi([0 1],1000,1));
    tx{i} = resample(tx{i},OSRs(i),1)/OSRs(i);
    tx{i} = hCarrierAggregationModulate(tx{i},SR,F_c(i)*1e6);
end

% Superpose the component carriers
waveform = tx{1};
for i = 2:numCC
    waveform = waveform + tx{i};
end

%%
specPlot = hCarrierAggregationPlotSpectrum(waveform,SR,...
    'Power Spectrum of Carrier Aggregation 4x20MHz R9',{'Signal spectrum'});