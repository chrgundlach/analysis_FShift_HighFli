%% parameters
clearvars
% F.Pathlocal             = 'G:\work\data\SSVEP_FShift_Probabil\';
F.Pathlocal             = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2025_FShiftHiFli\';


F.PathInEEG             = fullfile(F.Pathlocal, 'eeg\epoch\');
F.PathInBehavior        = fullfile(F.Pathlocal, 'behavior\');
F.PathInSCADS           = fullfile(F.Pathlocal, 'eeg\SCADS\');
F.PathOut               = fullfile(F.Pathlocal, 'eeg\fft\'); % with FWHM 0.5
F.subjects              = arrayfun(@(x) sprintf('%02.0f',x),1:70,'UniformOutput',false)';
% F.sub2use               = [6:13 15:18];%:53;
% F.sub2use               = [1:8];%
F.sub2use               = [21 22];%

% changed experiment from participant 22 onwards (stimuli isoluminant to background

F.trigger               =   {[10 ]; ... %RDK1  attended conventional flicker
                             [20 ]; ... %RDK2  attended conventional flicker
                             [30 ]; ... %RDK1 attended HighFreqFlick I colored inlet
                             [40 ]; ... %RDK2 attended HighFreqFlick I colored inlet
                             [50 ]; ... %RDK1 attended HighFreqFlick II white inlet
                             [60 ]}; ... %RDK2 attended HighFreqFlick II white inlet

F.EEGChans              = 64;

F.Trials2Consider       = [0 1]; % trials to consider, i.e.  [0 1] all; [0.5 1] second half
F.TotTrialNum           = 360; % total number of trials [without events]

F.FFT_epoch             = [-1 2];
F.CSD_flag              = 1; % 0 = no; 1 = yes

F.FFT_timewins          = {[-1 0]; [0.5 1.5]; [1 2]}; % for v1
% F.FFT_timewins          = {[-0.5 0]; [0.25 0.75]; [0.75 1.25];[1 1.5];[1.25 1.75]}; % for v2
% F.FFT_timewins          = {[-0.75 0]; [0.25 1]; [1 1.75]}; % for v3
F.FFT_freqres           = 16384;



%% start processing
%% loop across subjects
for i_sub = 1:numel(F.sub2use)
    %% load files
    % EEG
    EEG = pop_loadset('filename',sprintf('VP%s_e.set',F.subjects{F.sub2use(i_sub)}),'filepath',F.PathInEEG);
    Preprocessing = load(fullfile(F.PathInSCADS,sprintf('VP%s_Preprocess_summary.mat',F.subjects{F.sub2use(i_sub)})));
    % pop_eegplot(EEG,1,1,1)
    % behavior (loads latest file)
    t.files = dir(fullfile(F.PathInBehavior,sprintf('VP%s_timing*.mat',F.subjects{F.sub2use(i_sub)})));
    [t.val t.idx ]=max([t.files.datenum]);
    behavior = load(fullfile(F.PathInBehavior,t.files(t.idx).name));
    
    % select subset of trials
%     t.alltrials = 1:F.TotTrialNum; % vector with all trials
%     t.alltrials(~(Preprocessing.PreProc.trial_blink & Preprocessing.PreProc.trial_eyemov & Preprocessing.PreProc.trial_SCADS))=[];
%     t.trialnumremapped = linspace(0,1,F.TotTrialNum); % map num of trials to [0 1]
%     t.idx = find(ismember(... % find trials that fall in the range of F.Trials2Consider
%         t.alltrials, find(t.trialnumremapped >= F.Trials2Consider(1) & t.trialnumremapped <= F.Trials2Consider(2))));
%     EEG = pop_select(EEG,'trial',t.idx); % only select those trials    
    
    % select certain epochs
    t.idx = [EEG.event(ismember([EEG.event.type],unique(cell2mat(F.trigger)))).epoch];
    EEG = pop_select(EEG,'trial',t.idx);
    
    % select certain time window
    EEG = pop_select(EEG, 'time', F.FFT_epoch);
    
    %% do csd transform
    if F.CSD_flag==1
        if  i_sub == 1 % calculate CSD matrix
            try
                CSD.chanmat=ExtractMontage('C:\Dropboxdata\Dropbox\work\matlab\software\toolboxes\CSD\resource\10-5-System_Mastoids_EGI129.csd',{EEG.chanlocs.labels}');
            catch
                CSD.chanmat=ExtractMontage('C:\Users\EEG\Documents\MATLAB\christopher\general_functions\CSD\resource\10-5-System_Mastoids_EGI129.csd',{EEG.chanlocs.labels}');
            end
            [CSD.G,CSD.H] = GetGH(CSD.chanmat);
        end
        fprintf(1,'\n###\ncalculating CSD transform\n###\n')
        for i_tr = 1:EEG.trials
            % csd of raw data
            EEG.data(:,:,i_tr)= CSDTransform(EEG.data(:,:,i_tr), CSD.G, CSD.H);
        end
    end
    
    
    %% some bookkeeping
    FFT.electrodes = EEG.chanlocs;
    FFT.RDK = behavior.RDK;
    
    FFT.params.srate = EEG.srate;
    
    %filter
    FFT.params.filter = {[] []};
    
    
    FFT.all_trialnum = EEG.trials;
    FFT.con_trigger = F.trigger;
    FFT.con_trialnum = cellfun(@(x) sum(ismember([EEG.event.type],x)), F.trigger);
    
    % index trials
    t.trialindex = cellfun(@(x) [EEG.event([ismember([EEG.event.type],x)]).epoch], F.trigger,'UniformOutput',false);
    
        
    %% save FFT and coherence transforms
    FFT.FFT.res = F.FFT_freqres;
    FFT.FFT.timewin = F.FFT_timewins;
    FFT.FFT.amp_data_ind=nan(FFT.FFT.res, F.EEGChans,numel(F.trigger),numel(FFT.FFT.timewin));
    FFT.FFT.amp_data_evo= FFT.FFT.amp_data_ind;
    FFT.FFT.mscohere_data_evo= FFT.FFT.amp_data_ind;
    FFT.FFT.PLV_data_evo= FFT.FFT.amp_data_ind;
    fprintf('FFT: %1.0f windows for %1.0f electrodes:\n', numel(FFT.FFT.timewin), EEG.nbchan)
    for i_win = 1:numel(FFT.FFT.timewin)
        % first: select relevant eeg data
        EEGt = pop_select(EEG,'time',FFT.FFT.timewin{i_win});
        % surrogate mean data for evoked spectra
        EEGt_m = pop_select(EEGt,'trial',1);
        % detrend
        EEGt = eegF_Detrend(EEGt,[]);
        
        % extract complex loop across channels
        fprintf('FFT for channel:  ')
        t.FFTdata = nan(FFT.FFT.res, F.EEGChans,EEGt.trials);
        for i_el = 1:EEGt.nbchan
            fprintf('\b\b%02.0f',i_el)
            t.FFTdata(:,i_el,:) = squeeze(fft(EEGt.data(i_el,:,:),FFT.FFT.res,2));
        end
        FFT.FFT.freqs = ((0:size(t.FFTdata,1)-1)/size(t.FFTdata,1)) * EEGt.srate;
        fprintf('\nCoherence+Amplitude+PLV for conditions:')
        % loop across conditions
        for i_con = 1:numel(FFT.con_trigger)
            fprintf('\b\b%02.0f',i_con)
            % amplitude values from FFT
            % induced amplitude
            FFT.FFT.amp_data_ind(:,:,i_con,i_win)=mean( ...
                abs(t.FFTdata(:,:,t.trialindex{i_con}))*2/size(EEGt.data,2) ...
                ,3);
            % figure; plot(FFT.FFT.freqs, FFT.FFT.amp_data_ind(:,29,i_con,i_win))
            
            % evoked amplitude by averaging trial data first 
            EEGt_m.data = detrend(mean(EEGt.data(:,:,t.trialindex{i_con}),3)')';
            % evoked amplitude from trial-averaged detrended data
            FFT.FFT.amp_data_evo(:,:,i_con,i_win) = squeeze(abs(fft(EEGt_m.data,FFT.FFT.res,2))*2/size(EEGt.data,2))';
            % figure; plot(FFT.FFT.freqs, FFT.FFT.data_evo(:,29,i_con,i_win))

            % PLV based on complex fft
            FFT.FFT.PLV_data_evo(:,:,i_con,i_win) = ...
                abs(mean(t.FFTdata(:,:,t.trialindex{i_con}) ./ abs(t.FFTdata(:,:,t.trialindex{i_con})), 3));
            % figure; plot(FFT.FFT.freqs, FFT.FFT.PLV_data_evo(:,29,i_con,i_win))

                       
            % magnitude squared coherence based on complex fft
            % or rather Ratio of Evoked Power to Total Power, often referred to as Coherent Power Fraction or Inter-Trial Linear Coherence (Squared)
            %coh(chans,t) = (  abs(sum(mX(:,t).*mY(:,t).*exp(1).^(1i*phase_diff(:,t)))/n)  .^2)  /  (sum(  (mX(:,t).^2).*(mY(:,t).^2)  )/n);
            FFT.FFT.mscohere_data_evo(:,:,i_con,i_win) = ...
                (abs( ...
                sum(t.FFTdata(:,:,t.trialindex{i_con}),3) ... 
                /size(t.FFTdata(:,:,t.trialindex{i_con}),3) ...
                ).^2) ./ ...
                (sum(abs(t.FFTdata(:,:,t.trialindex{i_con})).^2,3)/size(t.FFTdata(:,:,t.trialindex{i_con}),3));
            % figure; plot(FFT.FFT.freqs, FFT.FFT.mscohere_data_evo(:,29,i_con,i_win))
           
            
            
        end
        fprintf('\n')
    end

    
    %% save
    FFT.savetime=datestr(now);
    if ~exist(F.PathOut); mkdir(F.PathOut); end
    fprintf(1,'|| saving file ||  %s\\VP%s_FFT.mat ||\n', ...
        F.PathOut,F.subjects{F.sub2use(i_sub)})
    save(fullfile(F.PathOut, sprintf('VP%s_FFT.mat',F.subjects{F.sub2use(i_sub)})), 'FFT')
    
end