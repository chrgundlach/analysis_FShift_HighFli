%% script to analyze behavioral data
clearvars

p.path =                'N:\AllgPsy\experimental_data\2025_FShift_Prime1of2\behavior\';
p.path =                '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2025_FShiftHiFli\Behavior\';
p.subs=                 arrayfun(@(x) sprintf('%02.0f',x),1:70,'UniformOutput',false)';
% p.subs2use=             [1:13 15:21];% 
% changed experiment from participant 22 onwards (stimuli isoluminant to background
p.subs2use=             [1:20];%

% p.subs2use=             [1:10];% %VP 15 no data
p.responsewin =         [0.2 1]; % according to p.targ_respwin from run_posnalpha
% p.responsewin =         [0.9 1.3]; % according to p.targ_respwin from run_posnalpha
p.resprecalcflag =      true; %do recalculation

p.eventtype =           {'target';'distractor'};
% p.eventtype2 =          {'target_primed';'target_nonprimed';'distractor'};


p.flickertype =        {'conventional';'conventional';'highfreq_col';'highfreq_col';'highfreq_white';'highfreq_white'};

p.colnames2change =     {'blue', 'green'};



%% actual calculation
for i_sub = 1:numel(p.subs2use)
    % load data
    temp.files = dir(sprintf('%sVP%s_timing*.mat',p.path,p.subs{p.subs2use(i_sub)}));
    data_in.resp.experiment = repmat({[nan]},1,12);
    data_in.button_presses.experiment = repmat({[nan]},1,12);
    for i_fi = 1:numel(temp.files)
        temp.data_in{i_fi}=load(sprintf('%s%s',p.path,temp.files(i_fi).name));
        % extract relevant data
        try data_in.conmat.experiment = temp.data_in{i_fi}.randmat.experiment;
            data_in.RDK = temp.data_in{i_fi}.RDK.RDK;
        end
        if any(strcmp(fieldnames(temp.data_in{i_fi}.resp),'experiment'))
            temp.index=cell2mat(cellfun(@(x) ~isempty(cell2mat({x})), temp.data_in{i_fi}.resp.experiment,'UniformOutput',false));
            data_in.resp.experiment(temp.index)=temp.data_in{i_fi}.resp.experiment(temp.index);
            data_in.button_presses.experiment(temp.index)=temp.data_in{i_fi}.button_presses.experiment(temp.index);
        end
    end
    % extract aprameter
    params(i_sub) = temp.data_in{1}.p;
    %% recalculate events
    if any(params(i_sub).targ_respwin-(p.responsewin*1000)~=0) & p.resprecalcflag
        % loop across blocks
        for i_block = 1:numel(data_in.resp.experiment)
            % index relevant response button
            t.buttonidx = data_in.button_presses.experiment{i_block}.keymap_ind == data_in.button_presses.experiment{i_block}.SPACE;
            % loop across trials
            for i_trial=1:numel(data_in.resp.experiment{i_block})
                % index trial data
                t.trial = data_in.resp.experiment{i_block}(i_trial);
                % index all button presses
                t.response = data_in.button_presses.experiment{i_block}.presses{i_trial}(:,t.buttonidx);
                % index time of button presses
                t.response_t = data_in.button_presses.experiment{i_block}.presses_t{i_trial}(:,t.buttonidx);
                % check for events?
                if any(~isnan(t.trial.eventtype))
                    t.idx = find(~isnan(t.trial.eventtype));
                    % loop across events
                    for i_eve = 1:numel(t.idx)
                        % index onset times/frames
                        t.onset_t = t.trial.event_onset_times(t.idx(i_eve));
                        t.onset_fr = t.trial.event_onset_frames(t.idx(i_eve));

                        % check whether any button press falls within response time window
                        t.respidx = find(t.response_t > t.onset_t+p.responsewin(1) & t.response_t < t.onset_t+p.responsewin(2));
                        % calculate RT
                        t.rt = (t.response_t(t.respidx)-t.onset_t)*1000;

                        % switch between distractor and event
                        if t.trial.eventtype(i_eve) == 1 % target
                            % response?
                            if ~isempty(t.rt)
                                data_in.resp.experiment{i_block}(i_trial).event_response_type{i_eve}='hit';
                                data_in.resp.experiment{i_block}(i_trial).event_response_RT(i_eve)=t.rt(1);
                            else
                                data_in.resp.experiment{i_block}(i_trial).event_response_type{i_eve}='miss';
                                data_in.resp.experiment{i_block}(i_trial).event_response_RT(i_eve)=nan;
                            end
                        else % distractor
                            if ~isempty(t.rt)
                                data_in.resp.experiment{i_block}(i_trial).event_response_type{i_eve}='FA_proper';
                                data_in.resp.experiment{i_block}(i_trial).event_response_RT(i_eve)=t.rt(1);
                            else
                                data_in.resp.experiment{i_block}(i_trial).event_response_type{i_eve}='CR';
                                data_in.resp.experiment{i_block}(i_trial).event_response_RT(i_eve)=nan;
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% setup response matrices
    data_all_struct = horzcat(data_in.resp.experiment{:});
    data_all_table = struct2table(data_all_struct);
    clear data_events
    
    % extract single event data
    t.evidx = ~isnan([data_all_struct.eventtype]);
    [t.data_eventtype]= deal(repmat({''},size([data_all_struct.eventtype])));
    t.mat = [data_all_struct.eventtype];
    t.data_eventtype(t.evidx) = arrayfun(@(x) p.eventtype{x}, t.mat(~isnan(t.mat)),'UniformOutput',false);
    t.pre_cuetimes = [data_all_struct.pre_cue_times];
    t.data_eventonsettimes = [data_all_struct.event_onset_times] - repmat(t.pre_cuetimes,2,1);
    t.data_trialnum = repmat([data_all_struct.trialnumber],2,1); % trialnumber
    t.data_blocknum = repmat([data_all_struct.blocknumber],2,1); % block number
    t.data_condition = repmat([data_all_struct.condition],2,1); % condition
    t.data_flickertype = arrayfun(@(x) p.flickertype(x),t.data_condition);
    t.data_cuedRDK = repmat([data_all_struct.cue],2,1); % condition
    t.data_eventRDK = [data_all_struct.eventRDK]; % RDK of event
    t.data_eventcolorname = repmat({''},size(t.data_eventRDK));
    for i_rdk = 1:numel(data_in.RDK)
        if strcmp(p.colnames2change(1,1),data_in.RDK(i_rdk).colnames)
            t.data_eventcolorname(t.data_eventRDK==i_rdk)=p.colnames2change(1,2);
        else
            t.data_eventcolorname(t.data_eventRDK==i_rdk)=deal({data_in.RDK(i_rdk).colnames});
        end
    end
    t.data_eventfreq = nan(size(t.data_eventRDK));
    % index differerent event frequencies
    t.idx = t.data_condition <=2; % conventional flicker
    t.data_eventfreq(t.evidx&t.idx)=arrayfun(@(x) data_in.RDK(x).freq, t.data_eventRDK(t.evidx&t.idx));
    t.idx = t.data_condition >2; % inlet flicker
    t.data_eventfreq(t.evidx&t.idx)=arrayfun(@(x) data_in.RDK(x).in_freq, t.data_eventRDK(t.evidx&t.idx));
    t.data_eventnum = repmat([1;2],1,size(t.evidx,2)).*t.evidx;

    t.data_response = repmat({''},size(t.data_eventRDK));
    t.data_response(~isnan(t.data_eventRDK)) = [data_all_struct.event_response_type];
    t.data_response_RT = nan(size(t.data_eventRDK));
    t.data_response_RT(~isnan(t.data_eventRDK))  = [data_all_struct.event_response_RT];
    t.data_subject = repmat(p.subs2use(i_sub),size(t.data_eventtype));
    %t.data_irr_color = repmat({data_in.RDK(5).colnames},size(t.data_eventtype));
    %t.data_irr_freq = repmat(data_in.RDK(3).freq,size(t.data_eventtype));
    

   
    
    % write into one file
    data_events = cell2struct([num2cell([t.data_subject(t.evidx) t.data_trialnum(t.evidx) t.data_blocknum(t.evidx) ...
        t.data_condition(t.evidx)]) t.data_flickertype(t.evidx) t.data_eventtype(t.evidx) ...
        num2cell(t.data_eventRDK(t.evidx)) t.data_eventcolorname(t.evidx) ...
        num2cell([t.data_eventfreq(t.evidx) t.data_eventnum(t.evidx) t.data_eventonsettimes(t.evidx)]) ...
        t.data_response(t.evidx) num2cell(t.data_response_RT(t.evidx))]',...
        {'participant','trialnumber','blocknumber', ...
        'condition','flickertype','eventtype', ...
        'RDKnumber','eventcolor', ...
        'eventfrequency','eventnumber','eventonsettimes',...
        'response','RT'});
    
    if i_sub == 1
        response_events = data_events;
    else
        response_events = vertcat(response_events,data_events);
    end
    
    %% check for false alarms?
    response_FA(i_sub).subject = p.subs2use(i_sub);
    response_FA(i_sub).FA_proper = sum(strcmpi({data_events.response},'FA_proper'));
    response_FA(i_sub).CR = sum(strcmpi({data_events.response},'CR'));
    response_FA(i_sub).FA_proper_rate = sum(strcmpi({data_events.response},'FA_proper'))./ ( ...
        sum(strcmpi({data_events.response},'FA_proper')) + sum(strcmpi({data_events.response},'CR')));
    response_FA(i_sub).FA = ...
        sum(cellfun(@(x) sum(cellfun(@(y) sum(y(:,8)~=0),x.presses_t)),data_in.button_presses.experiment))-...
        sum(strcmpi({data_events.response},'hit')|strcmpi({data_events.response},'FA_proper'));
    
end

%% interim statistics
response_events_table = struct2table(response_events);
t.dat = grpstats(response_events_table, ["participant","response"],'numel');
t.dat2 = grpstats(response_events_table, ["participant","response"],'mean',"DataVars",["RT"]);

t.idx = strcmp(t.dat.response,'FA_proper');
[t.dat.participant(1:4:end) t.dat.GroupCount(strcmp(t.dat.response,'FA_proper')) t.dat.GroupCount(strcmp(t.dat.response,'CR')) ...
    100*[t.dat.GroupCount(strcmp(t.dat.response,'FA_proper'))./(t.dat.GroupCount(strcmp(t.dat.response,'CR')) + ...
    t.dat.GroupCount(strcmp(t.dat.response,'FA_proper')))]]
mean(100*[t.dat.GroupCount(strcmp(t.dat.response,'FA_proper'))./(t.dat.GroupCount(strcmp(t.dat.response,'CR')) + ...
    t.dat.GroupCount(strcmp(t.dat.response,'FA_proper')))])

t.dat.GroupCount(strcmp(t.dat.response,'hit'))./(t.dat.GroupCount(strcmp(t.dat.response,'hit'))+t.dat.GroupCount(strcmp(t.dat.response,'miss')))

mean(t.dat2.mean_RT(strcmp(t.dat.response,'hit')))

%% summary statistics/descriptives
response_events_table = struct2table(response_events);
response_FA_table = struct2table(response_FA);
% head(response_events_table)

% general summary mean RT and hit rate
% export data
p.path = 'C:\Dropboxdata\Dropbox\work\R-statistics\experiments\ssvep_fshifthighfreqflick\data_in';
% p.path = 'C:\Users\EEG\Documents\R\Christopher\analysis_R_ssvep_fshift_perirr\data_in';
p.filename1 = 'behavior_events.csv';
p.filename2 = 'behavior_FAs.csv';
writetable(response_events_table, fullfile(p.path,p.filename1))
writetable(response_FA_table, fullfile(p.path,p.filename2))

