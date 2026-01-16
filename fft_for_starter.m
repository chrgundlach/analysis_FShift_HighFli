%% script to start with fft analysis

clearvars

%% parameters
% F.Pathlocal             = 'G:\work\data\SSVEP_FShift_Probabil\';
F.Pathlocal             = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2025_FShiftHiFli\';


F.subjects              = arrayfun(@(x) sprintf('%02.0f',x),1:70,'UniformOutput',false)';
F.PathInEEG             = fullfile(F.Pathlocal, 'eeg\epoch\');
F.sub2use               = [1];%

% changed experiment from participant 22 onwards (stimuli isoluminant to background

F.trigger               =   {[10 ]; ... %RDK1/3  attended conventional flicker
                             [20 ]; ... %RDK2/4  attended conventional flicker
                             [30 ]; ... %RDK1 attended HighFreqFlick I colored inlet
                             [40 ]; ... %RDK2 attended HighFreqFlick I colored inlet
                             [50 ]; ... %RDK1 attended HighFreqFlick II white inlet
                             [60 ]}; ... %RDK2 attended HighFreqFlick II white inlet

F.con_flickertype = [1 1 2 2 3 3];
F.conlabel_flicker = {'conventional';'conventional';'highfreq_col';'highfreq_col';'highfreq_white';'highfreq_white'};
F.con_flickerfreq = [20 23 14 17; 20 23 14 17; 68 71 62 65 ; 68 71 62 65; 68 71 62 65; 68 71 62 65];

i_sub = 1
%% load files
% EEG
EEG = pop_loadset('filename',sprintf('VP%s_e.set',F.subjects{F.sub2use(i_sub)}),'filepath',F.PathInEEG);

%% select the relevant trials
t.con2plot = [5 6];
t.idx = [EEG.event(ismember([EEG.event.type],unique(cell2mat(F.trigger(t.con2plot))))).epoch];
EEG_trialsofinterest = pop_select(EEG,'trial',t.idx); % select trials if interest

% time window of interst
t.time2plot = [-1 2]; % time in s
EEG_trialsandtimeofinterest = pop_select(EEG_trialsofinterest,'time',t.time2plot);
pop_eegplot(EEG_trialsandtimeofinterest,1,1,1)

% plot Oz data
t.electrode = 29;
figure; plot(EEG_trialsandtimeofinterest.times,squeeze(EEG_trialsandtimeofinterest.data(t.electrode,:,:)))
hold on
plot(EEG_trialsandtimeofinterest.times, ...
    squeeze(mean(EEG_trialsandtimeofinterest.data(t.electrode,:,:),3)),'k','LineWidth',2)
xlabel('time in ms')
ylabel('amplitude in \muV')
title('the best Oz signal ever | conventional flicker | single trials + average')

% fft or what?
% Perform FFT on the Oz data
% phaselocking!
Fs = EEG_trialsandtimeofinterest.srate; % Sampling frequency
dataelec = squeeze(mean(EEG_trialsandtimeofinterest.data(t.electrode,:,:),3)); % Extract Oz data
n = size(dataelec, 2); % Number of samples
f = (0:n-1)*(Fs/n); % Frequency range
Y = fft(dataelec, [], 2); % Compute FFT along the time dimension
Amp = abs(Y)*2/n; % amplitude of frequencies mu/V

% plot this
figure;
plot(f, Amp)
ylabel('amplitude in \muV')
xlabel('frequency in Hz')
title('plot this FFT of Oz')
xlim([0 70])

%% we want to plot a topography
% for which frequencies?
t.f2plot = [68 71]; % central high frequencies

% fft for all electrodes
Fs = EEG_trialsandtimeofinterest.srate; % Sampling frequency
dataelec = squeeze(mean(EEG_trialsandtimeofinterest.data(:,:,:),3)); % Extract Oz data
n = size(dataelec, 2); % Number of samples
f = (0:n-1)*(Fs/n); % Frequency range
Y = fft(dataelec, [], 2); % Compute FFT along the time dimension
Amp = abs(Y)*2/n; % amplitude of frequencies mu/V


% Compute the topographic map for the specified frequencies
for f_idx = 1:length(t.f2plot)
    freq = t.f2plot(f_idx);
    % Find the index of the frequency in the FFT result
    [~, freq_idx] = min(abs(f - freq));
    topography(f_idx, :) = Amp( :,freq_idx); % Store amplitude for the frequency
end

figure;
topoplot(mean(topography,1),EEG.chanlocs,'maplimits',[0 max(mean(topography,1))]);
colorbar
% colormap("parula")