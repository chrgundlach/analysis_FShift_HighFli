%% plot FFT images
clearvars
F.PathInEEG             = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2025_FShiftHiFli\EEG\fft'; % with FWHM 0.5

F.Subs                  = arrayfun(@(x) sprintf('%02.0f',x),1:70,'UniformOutput',false)';
% F.Subs2use              = [1:13 15:21];
% changed experiment from participant 22 onwards (stimuli isoluminant to
% background and used other frequencies
% participant 42 has lower trial number
F.Subs2use              = [1:8]; 
                        
F.fft.baseline          = [-500 -250];

F.SSVEP_Freqs           = [20 23 14 17]; 
F.RDK_pos               = [0 0 -255 -255];
F.RDK_pos_label         = {'center';'center'; 'left';'left'};



F.conlabel_att = {'att RDK1';'att RDK2'; 'att RDK1';'att RDK2';'att RDK1';'att RDK2'};
F.conRDKattended = logical([1 0 1 0; 0 1 0 1;1 0 1 0; 0 1 0 1;1 0 1 0; 0 1 0 1]);
F.conRDKattended_label = repmat({'attended'},size(F.conRDKattended));
F.conRDKattended_label(F.conRDKattended==0) = {'not attended'};
F.con_flickertype = [1 1 2 2 3 3];
F.conlabel_flicker = {'conventional';'conventional';'highfreq_col';'highfreq_col';'highfreq_white';'highfreq_white'};
F.con_flickerfreq = [20 23 14 17; 20 23 14 17; 68 71 62 65 ; 68 71 62 65; 68 71 62 65; 68 71 62 65];
F.con_RDKpresented = [1 2 3 4; 1 2 3 4; 5 6 7 8; 5 6 7 8; 9 10 11 12; 9 10 11 12];

pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];
%% load data
for i_sub = 1:numel(F.Subs2use)
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%s_FFTdat.mat ||\n',...
        i_sub,numel(F.Subs2use),F.PathInEEG,F.Subs{F.Subs2use(i_sub)})
    
    temp.fft = open(fullfile(F.PathInEEG,sprintf('VP%s_fft.mat',F.Subs{F.Subs2use(i_sub)})));
    
    % convert to single
    temp.fft.FFT.FFT.amp_data_ind = single(temp.fft.FFT.FFT.amp_data_ind);
    temp.fft.FFT.FFT.amp_data_evo = single(temp.fft.FFT.FFT.amp_data_evo);
    temp.fft.FFT.FFT.mscohere_data_evo = single(temp.fft.FFT.FFT.mscohere_data_evo);
    temp.fft.FFT.FFT.PLV_data_evo = single(temp.fft.FFT.FFT.PLV_data_evo);
    
    
    % preallocate memory
    if i_sub == 1
        FFTdat.electrodes = temp.fft.FFT.electrodes;
        FFTdat.con_trialnum = temp.fft.FFT.con_trialnum;
        FFTdat.srate = temp.fft.FFT.params.srate/2;
        FFTdat.amp_data_evo = nan([size(temp.fft.FFT.FFT.amp_data_evo),numel(F.Subs2use)]);
        FFTdat.amp_data_ind = nan([size(temp.fft.FFT.FFT.amp_data_ind),numel(F.Subs2use)]);
        FFTdat.mscohere_data_evo = nan([size(temp.fft.FFT.FFT.mscohere_data_evo),numel(F.Subs2use)]);
        FFTdat.PLV_data_evo = nan([size(temp.fft.FFT.FFT.PLV_data_evo),numel(F.Subs2use)]);
        FFTdat.ffttimewin = temp.fft.FFT.FFT.timewin;
        FFTdat.fftfreqs = temp.fft.FFT.FFT.freqs;

    end
    
    % assign data
    FFTdat.amp_data_evo(:,:,:,:,i_sub) = temp.fft.FFT.FFT.amp_data_evo; % evoked data
    FFTdat.amp_data_ind(:,:,:,:,i_sub) = temp.fft.FFT.FFT.amp_data_ind; % induced data
    FFTdat.mscohere_data_evo(:,:,:,:,i_sub) = temp.fft.FFT.FFT.mscohere_data_evo;
    FFTdat.PLV_data_evo(:,:,:,:,i_sub) = temp.fft.FFT.FFT.PLV_data_evo;
    FFTdat.RDK(i_sub) = temp.fft.FFT.RDK;

      
    clear temp    
end

pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];



%% plot electrode head
% pl.elec2plot = {'C3';'CP3'};
pl.elec2plot = {'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'P8';'P10';'PO4';'PO8';'O2';'I2'};
% pl.elec2plot = {'P6';'P8';'P10';'PO4';'PO8';'O2';'I2';'POz';'Oz';'Iz';'O1'};
% pl.elec2plot = {'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'O2'};
% pl.elec2plot = {'POz';'O1';'Oz';'O2';'Iz'};
% pl.elec2plot = {'P8';'P10';'PO8';'PO3';'POz';'Oz';'O1'};
% pl.elec2plot = {'P7';'P9';'PO7';'PO2';'POz';'Oz';'O2'};


figure;
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({FFTdat.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
% topoplot(find(pl.elec2plot_i),FFTdat.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on');
topoplot([],FFTdat.electrodes(1:64),'whitebk','on','style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i),'o','r',8});

%% plot grand mean FFT data | spectra | for subset of conditions

% large center as in tango | periphery: central and lateral 
% pl.elec2plot = {{'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}, 'center'};
% pl.elec2plot = {{'P7';'P5';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P8';'P6';'PO4';'PO8';'P10';'O2';'I2'}, 'center'};
% pl.elec2plot = {{'POz';'Oz';'Iz';'P8';'P6';'PO4';'PO8';'P10';'O2';'I2'}, 'left'}; % topo peak in right hemisphere?
% pl.elec2plot = { {'PO3';'POz';'O1';'Oz';'PO8';'P8';'P10'}, 'left'}; % topo peak in right hemisphere?


% pl.elec2plot = {{'PO3';'POz';'PO4';'O1';'Oz';'O2'}, 'center'}; % topo peak in right hemisphere?
% pl.elec2plot = {{'PO3';'POz';'PO4';'O1';'Oz';'O2';'I1';'Iz';'I2'}, 'center'}; % topo peak in right hemisphere?
pl.elec2plot = {{'PO3';'POz';'PO4';'O1';'Oz';'O2';'I1';'Iz';'I2'}, 'left'}; % topo peak in right hemisphere?

% pl.elec2plot = {{'P5';'P3';'PO7';'PO3';'P4';'P6';'PO4';'PO8';'O1';'Oz';'O2';'I1';'Iz';'I2'}, 'center SSVEP conventional'
%     {'PO3';'POz';'O1';'Oz';'PO8';'P8';'P10'}, 'left SSVEP conventional';
%     {'PO3';'POz';'PO4';'O1';'Oz';'O2'}, 'SSVEP high freq'};


pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({FFTdat.electrodes.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);

% pl.time2plot = [1];
pl.time2plot = [1:3];
pl.sub2plot = 1:numel(F.Subs2use);

% pl.con2plot = [1 2];
% pl.con2plot = [3 4];
pl.con2plot = [5 6];

% extract data
pl.data_ind = squeeze(mean(FFTdat.amp_data_ind(:,pl.elec2plot_i{1},pl.con2plot,pl.time2plot,pl.sub2plot),[2,3,4]));
pl.data_evo = squeeze(mean(FFTdat.amp_data_evo(:,pl.elec2plot_i{1},pl.con2plot,pl.time2plot,pl.sub2plot),[2,3,4]));
pl.data_mscohere = squeeze(mean(FFTdat.mscohere_data_evo(:,pl.elec2plot_i{1},pl.con2plot,pl.time2plot,pl.sub2plot),[2,3,4]));
pl.data_PLV = squeeze(mean(FFTdat.PLV_data_evo(:,pl.elec2plot_i{1},pl.con2plot,pl.time2plot,pl.sub2plot),[2,3,4]));

% FOI determined by cons and position
pl.FOI = F.con_flickerfreq(pl.con2plot(1),strcmp(F.RDK_pos_label,pl.elec2plot{2}));

% plotting
figure;
set(gcf,'Position',[100 100 1100 700],'PaperPositionMode','auto')
subplot(2,2,1)
plot(FFTdat.fftfreqs,pl.data_ind,'Color',[0.5 0.5 0.5],'LineWidth',1)
hold on;
plot(FFTdat.fftfreqs,mean(pl.data_ind,2),'Color','k','LineWidth',2)
xlim([0 75])
ylim(get(gca,'ylim'))
xlabel('frequency in Hz')
ylabel('amplitude in \muV/m²')
title(sprintf('induced GrandMean FFT spectra | %s\nN = %1.0f | FOI = %1.1f %1.1f  Hz', ...
    F.conlabel_flicker{pl.con2plot(1)},numel(pl.sub2plot),pl.FOI),'Interpreter','none')
vline(pl.FOI,'k:')

subplot(2,2,3)
hold on;
plot(FFTdat.fftfreqs,pl.data_evo,'Color',[0.5 0.5 0.5],'LineWidth',1)
plot(FFTdat.fftfreqs,mean(pl.data_evo,2),'Color','k','LineWidth',2)
xlim([0 75])
ylim(get(gca,'ylim'))
xlabel('frequency in Hz')
ylabel('amplitude in \muV/m²')
title(sprintf('evoked GrandMean FFT spectra | %s\nN = %1.0f | FOI = %1.1f %1.1f Hz', ...
     F.conlabel_flicker{pl.con2plot(1)},numel(pl.sub2plot), pl.FOI),'Interpreter','none')
vline(pl.FOI,'k:')
box on

subplot(2,2,2)
hold on;
plot(FFTdat.fftfreqs,pl.data_PLV,'Color',[0.5 0.5 0.5],'LineWidth',1)
plot(FFTdat.fftfreqs,mean(pl.data_PLV,2),'Color','k','LineWidth',2)
xlim([0 75])
ylim(get(gca,'ylim'))
xlabel('frequency in Hz')
ylabel('PLV')
title(sprintf('evoked GrandMean PLV spectra | %s\nN = %1.0f | FOI = %1.1f %1.1f Hz', ...
     F.conlabel_flicker{pl.con2plot(1)},numel(pl.sub2plot), pl.FOI),'Interpreter','none')
vline(pl.FOI,'k:')
box on

subplot(2,2,4)
hold on;
plot(FFTdat.fftfreqs,pl.data_mscohere,'Color',[0.5 0.5 0.5],'LineWidth',1)
plot(FFTdat.fftfreqs,mean(pl.data_mscohere,2),'Color','k','LineWidth',2)
xlim([0 75])
ylim(get(gca,'ylim'))
xlabel('frequency in Hz')
ylabel('PLV')
title(sprintf('evoked GrandMean mscohere spectra | %s\nN = %1.0f | FOI = %1.1f %1.1f Hz', ...
     F.conlabel_flicker{pl.con2plot(1)},numel(pl.sub2plot), pl.FOI),'Interpreter','none')
vline(pl.FOI,'k:')
box on


% draw topography with electrode positions
h.a1 = axes('position',[0 0.45 0.15 0.15],'Visible','off');
topoplot(find(any(cell2mat(pl.elec2plot_i))),FFTdat(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(any(cell2mat(pl.elec2plot_i),1)),'o','r',4,1});



%% plot Grand Mean FFT data | topoplot for different frequencies
pl.time2plot = [1:3];
pl.time2plot = [1];
% pl.freq2plot = [1 0 0 0]==1;
pl.freq2plot = [0 1 0 0]==1;
% pl.freq2plot = [0 0 1 0]==1;
% pl.freq2plot = [0 0 0 1]==1;
pl.freqrange=[-0.1 0.1];
pl.sub2sel = 1:numel(F.Subs2use);

pl.con2plot = [1 2];
% pl.con2plot = [3 4];
% pl.con2plot = [5 6];

pl.FOI=F.con_flickerfreq(pl.con2plot(1),pl.freq2plot);



t.idx = dsearchn(FFTdat.fftfreqs',(pl.freqrange+pl.FOI)');
% pl.data_ind = squeeze(mean(FFTdat.fftdata_ind(t.idx(1):t.idx(2),:,:,pl.time2plot,pl.sub2sel),[1 3 4 ]));
pl.data_evo = squeeze(mean(FFTdat.amp_data_evo(t.idx(1):t.idx(2),:,pl.con2plot,pl.time2plot,pl.sub2sel),[1 3 4 ]));

% try to subtract general topo?
pl.FOI2 = [min(F.con_flickerfreq(pl.con2plot(1),:))-3 max(F.con_flickerfreq(pl.con2plot(1),:))+3];
t.idx2 = dsearchn(FFTdat.fftfreqs',[pl.freqrange+pl.FOI2(1) pl.freqrange+pl.FOI2(2)]');
pl.data_evo_ctr = squeeze(mean(FFTdat.amp_data_evo([t.idx2(1):t.idx2(2) t.idx2(3):t.idx2(4)],:,pl.con2plot,pl.time2plot,pl.sub2sel),[1 3 4 ]));
pl.data_evo = pl.data_evo-pl.data_evo_ctr;


figure;
set(gcf,'Position',[100 100 300 300],'PaperPositionMode','auto')

% % induced
% h.s(1) = subplot(1,2,1);
% pl.mdata = mean(pl.data_ind,2,'omitnan');
% pl.clim = [0 max(pl.mdata)];
% topoplot(pl.mdata, FFTdat.electrodes(1:64), ...
%     'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',fake_parula,...
%     'whitebk','on');
% title(sprintf('induced SSVEP amps for %1.1f +- [%1.1f  %1.1f] Hz | [%1.0f %1.0f]ms', ...
%     pl.FOI, pl.freqrange, min([FFTdat.ffttimewin{pl.time2plot}])*1000, max([FFTdat.ffttimewin{pl.time2plot}])*1000))
% colorbar

% h.s(2) = subplot(1,2,2);
pl.mdata = mean(pl.data_evo,2,'omitnan');
pl.clim = [0 max(pl.mdata)];
topoplot(pl.mdata, FFTdat.electrodes(1:64), ...
    'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',fake_parula,...
    'whitebk','on');
title(sprintf('SSVEP amps %1.1f +- [%1.1f  %1.1f] Hz\n [%1.0f %1.0f]ms | %s', ...
    pl.FOI, pl.freqrange, min([FFTdat.ffttimewin{pl.time2plot}])*1000, max([FFTdat.ffttimewin{pl.time2plot}])*1000, ...
    F.RDK_pos_label{pl.freq2plot}))
colorbar



%% plot Grand Mean FFT data | topoplot for all SSVEP frequencies
pl.time2plot = [1:3];
pl.time2plot = [1];
pl.con2plot = [1 2;3 4;5 6];
pl.pos2plot = {'center';'left'};
pl.freqrange=[-0.1 0.1];
pl.sub2sel = 1:numel(F.Subs2use);

pl.data_evo = []; pl.data_evo_subt = []; 
pl.mscohere_data_evo = []; pl.mscohere_data_evo_subt = []; 
pl.PLV_data_evo = []; pl.PLV_data_evo_subt = [];t.FOI = [];
for i_con = 1:size(pl.con2plot,1)
    for i_pos = 1:numel(pl.pos2plot)
        % frequencies of interest
        t.FOI(:,i_con,i_pos) = unique(F.con_flickerfreq(pl.con2plot(i_con,:),strcmp(F.RDK_pos_label,pl.pos2plot{i_pos})));
        t.FOI2 = [min(F.con_flickerfreq(pl.con2plot(i_con,:),:),[],"all")-3 max(F.con_flickerfreq(pl.con2plot(i_con,:),:),[],"all")+3];
        
        % index frequencies
        t.t = dsearchn(FFTdat.fftfreqs',reshape((pl.freqrange+t.FOI(:,i_con,i_pos))',1,[])');
        t.fidx = unique(cell2mat(arrayfun(@(x,y) x:y,t.t(1:2:end),t.t(2:2:end),'UniformOutput',false)));
        t.t2 = dsearchn(FFTdat.fftfreqs',reshape((pl.freqrange+t.FOI2)',1,[])');
        t.fidx2 = unique(cell2mat(arrayfun(@(x,y) x:y,t.t2(1:2:end),t.t2(2:2:end),'UniformOutput',false)));

        % extract data
        pl.data_evo(:,i_con,i_pos) = squeeze(mean(FFTdat.amp_data_evo(t.fidx,:,pl.con2plot(i_con,:),pl.time2plot,pl.sub2sel),[1 3 4 5]))';
        pl.data_evo_subt(:,i_con,i_pos) = pl.data_evo(:,i_con,i_pos) - ...
            squeeze(mean(FFTdat.amp_data_evo(t.fidx2,:,pl.con2plot(i_con,:),pl.time2plot,pl.sub2sel),[1 3 4 5]))';
        pl.mscohere_data_evo(:,i_con,i_pos) = squeeze(mean(FFTdat.mscohere_data_evo(t.fidx,:,pl.con2plot(i_con,:),pl.time2plot,pl.sub2sel),[1 3 4 5]))';
        pl.mscohere_data_evo_subt(:,i_con,i_pos) = pl.mscohere_data_evo(:,i_con,i_pos) - ...
            squeeze(mean(FFTdat.mscohere_data_evo(t.fidx2,:,pl.con2plot(i_con,:),pl.time2plot,pl.sub2sel),[1 3 4 5]))';
        pl.PLV_data_evo(:,i_con,i_pos) = squeeze(mean(FFTdat.PLV_data_evo(t.fidx,:,pl.con2plot(i_con,:),pl.time2plot,pl.sub2sel),[1 3 4 5]))';
        pl.PLV_data_evo_subt(:,i_con,i_pos) = pl.PLV_data_evo(:,i_con,i_pos) - ...
            squeeze(mean(FFTdat.PLV_data_evo(t.fidx2,:,pl.con2plot(i_con,:),pl.time2plot,pl.sub2sel),[1 3 4 5]))';
    end
end

% plot data amp evoked
figure;
set(gcf,'Position',[100 100 1100 400],'PaperPositionMode','auto')

h = []; pl.num = 1;
for i_pos = 1:numel(pl.pos2plot)
    for i_con = 1:size(pl.con2plot,1)
        h.s(i_con,i_pos) = subplot(numel(pl.pos2plot),size(pl.con2plot,1),pl.num);
        pl.clim = [0 max(pl.data_evo(:,i_con,:),[],"all")];
        topoplot(pl.data_evo(:,i_con,i_pos), FFTdat.electrodes(1:64), ...
            'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',fake_parula,...
            'whitebk','on');
        title(sprintf('evo SSVEP [%1.1f %1.1f] +- [%1.1f  %1.1f] Hz\n[%1.0f %1.0f]ms | %s | %s', ...
            t.FOI(:,i_con,i_pos), pl.freqrange, min([FFTdat.ffttimewin{pl.time2plot}])*1000, max([FFTdat.ffttimewin{pl.time2plot}])*1000, ...
            pl.pos2plot{i_pos}, F.conlabel_flicker{pl.con2plot(i_con,1)}), "FontSize",8, "Interpreter","none")
        colorbar
        pl.num = pl.num+1;
    end
end


% plot data amp evoked for which neighboring frequencies were subtracted
figure;
set(gcf,'Position',[100 100 1100 400],'PaperPositionMode','auto')

h = []; pl.num = 1;
for i_pos = 1:numel(pl.pos2plot)
    for i_con = 1:size(pl.con2plot,1)
        h.s(i_con,i_pos) = subplot(numel(pl.pos2plot),size(pl.con2plot,1),pl.num);
        pl.clim = [0 max(pl.data_evo_subt(:,i_con,:),[],"all")];
        topoplot(pl.data_evo_subt(:,i_con,i_pos), FFTdat.electrodes(1:64), ...
            'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',fake_parula,...
            'whitebk','on');
        title(sprintf('evo subtr SSVEP [%1.1f %1.1f] +- [%1.1f  %1.1f] Hz\n[%1.0f %1.0f]ms | %s | %s', ...
            t.FOI(:,i_con,i_pos), pl.freqrange, min([FFTdat.ffttimewin{pl.time2plot}])*1000, max([FFTdat.ffttimewin{pl.time2plot}])*1000, ...
            pl.pos2plot{i_pos}, F.conlabel_flicker{pl.con2plot(i_con,1)}), "FontSize",8, "Interpreter","none")
        colorbar
        pl.num = pl.num+1;
    end
end


% plot mscohere
figure;
set(gcf,'Position',[100 100 1100 400],'PaperPositionMode','auto')

h = []; pl.num = 1;
for i_pos = 1:numel(pl.pos2plot)
    for i_con = 1:size(pl.con2plot,1)
        h.s(i_con,i_pos) = subplot(numel(pl.pos2plot),size(pl.con2plot,1),pl.num);
        pl.clim = [0 max(pl.mscohere_data_evo(:,i_con,:),[],"all")];
        topoplot(pl.mscohere_data_evo(:,i_con,i_pos), FFTdat.electrodes(1:64), ...
            'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',fake_parula,...
            'whitebk','on');
        title(sprintf('mscohere SSVEP [%1.1f %1.1f] +- [%1.1f  %1.1f] Hz\n[%1.0f %1.0f]ms | %s | %s', ...
            t.FOI(:,i_con,i_pos), pl.freqrange, min([FFTdat.ffttimewin{pl.time2plot}])*1000, max([FFTdat.ffttimewin{pl.time2plot}])*1000, ...
            pl.pos2plot{i_pos}, F.conlabel_flicker{pl.con2plot(i_con,1)}), "FontSize",8, "Interpreter","none")
        colorbar
        pl.num = pl.num+1;
    end
end


% plot mscohere data for which neighboring frequencies were subtracted
figure;
set(gcf,'Position',[100 100 1100 400],'PaperPositionMode','auto')

h = []; pl.num = 1;
for i_pos = 1:numel(pl.pos2plot)
    for i_con = 1:size(pl.con2plot,1)
        h.s(i_con,i_pos) = subplot(numel(pl.pos2plot),size(pl.con2plot,1),pl.num);
        pl.clim = [0 max(pl.mscohere_data_evo_subt(:,i_con,:),[],"all")];
        topoplot(pl.mscohere_data_evo_subt(:,i_con,i_pos), FFTdat.electrodes(1:64), ...
            'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',fake_parula,...
            'whitebk','on');
        title(sprintf('mscohere subtr SSVEP [%1.1f %1.1f] +- [%1.1f  %1.1f] Hz\n[%1.0f %1.0f]ms | %s | %s', ...
            t.FOI(:,i_con,i_pos), pl.freqrange, min([FFTdat.ffttimewin{pl.time2plot}])*1000, max([FFTdat.ffttimewin{pl.time2plot}])*1000, ...
            pl.pos2plot{i_pos}, F.conlabel_flicker{pl.con2plot(i_con,1)}), "FontSize",8, "Interpreter","none")
        colorbar
        pl.num = pl.num+1;
    end
end


% plot PLV data
figure;
set(gcf,'Position',[100 100 1100 400],'PaperPositionMode','auto')

h = []; pl.num = 1;
for i_pos = 1:numel(pl.pos2plot)
    for i_con = 1:size(pl.con2plot,1)
        h.s(i_con,i_pos) = subplot(numel(pl.pos2plot),size(pl.con2plot,1),pl.num);
        pl.clim = [0 max(pl.PLV_data_evo(:,i_con,:),[],"all")];
        topoplot(pl.PLV_data_evo(:,i_con,i_pos), FFTdat.electrodes(1:64), ...
            'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',fake_parula,...
            'whitebk','on');
        title(sprintf('PLV SSVEP [%1.1f %1.1f] +- [%1.1f  %1.1f] Hz\n[%1.0f %1.0f]ms | %s | %s', ...
            t.FOI(:,i_con,i_pos), pl.freqrange, min([FFTdat.ffttimewin{pl.time2plot}])*1000, max([FFTdat.ffttimewin{pl.time2plot}])*1000, ...
            pl.pos2plot{i_pos}, F.conlabel_flicker{pl.con2plot(i_con,1)}), "FontSize",8, "Interpreter","none")
        colorbar
        pl.num = pl.num+1;
    end
end


% plot mscohere data for which neighboring frequencies were subtracted
figure;
set(gcf,'Position',[100 100 1100 400],'PaperPositionMode','auto')

h = []; pl.num = 1;
for i_pos = 1:numel(pl.pos2plot)
    for i_con = 1:size(pl.con2plot,1)
        h.s(i_con,i_pos) = subplot(numel(pl.pos2plot),size(pl.con2plot,1),pl.num);
        pl.clim = [0 max(pl.PLV_data_evo_subt(:,i_con,:),[],"all")];
        topoplot(pl.PLV_data_evo_subt(:,i_con,i_pos), FFTdat.electrodes(1:64), ...
            'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',fake_parula,...
            'whitebk','on');
        title(sprintf('PLV subtr SSVEP [%1.1f %1.1f] +- [%1.1f  %1.1f] Hz\n[%1.0f %1.0f]ms | %s | %s', ...
            t.FOI(:,i_con,i_pos), pl.freqrange, min([FFTdat.ffttimewin{pl.time2plot}])*1000, max([FFTdat.ffttimewin{pl.time2plot}])*1000, ...
            pl.pos2plot{i_pos}, F.conlabel_flicker{pl.con2plot(i_con,1)}), "FontSize",8, "Interpreter","none")
        colorbar
        pl.num = pl.num+1;
    end
end










%% extract amplitude values for FFT
% plotting parameters
pl.elec2plot = {{'P5';'P3';'PO7';'PO3';'P4';'P6';'PO4';'PO8';'O1';'Oz';'O2';'I1';'Iz';'I2'}, 'center SSVEP conventional'
    {'PO3';'POz';'O1';'Oz';'PO8';'P8';'P10'}, 'left SSVEP conventional';
    {'PO3';'POz';'PO4';'O1';'Oz';'O2'}, 'SSVEP high freq'};
pl.elec2plot = {{'P5';'P3';'PO7';'PO3';'P4';'P6';'PO4';'PO8';'O1';'Oz';'O2';'I1';'Iz';'I2'}, 'center SSVEP conventional'
    {'PO3';'POz';'O1';'Oz';'PO8';'P8';'P10'}, 'left SSVEP conventional';
    {'PO3';'POz';'PO4';'O1';'Oz';'O2';'I1';'Iz';'I2'}, 'SSVEP high freq'};

pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({FFTdat.electrodes.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);

% plot electrode positions for checking
figure; set(gcf,'Position',[100 100 700 200],'PaperPositionMode','auto')
for i_elecs = 1:size(pl.elec2plot ,1)
    subplot(1,size(pl.elec2plot ,1),i_elecs)
    topoplot([],FFTdat.electrodes(1:64),'whitebk','on','style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i{i_elecs}),'o','r',8});
    title(pl.elec2plot{i_elecs,2})
end

pl.elec4conXrdk = [1 1 2 2; 1 1 2 2; repmat(3,4,4)]; % which electrode cluster for which condition x RDK?

pl.time2plot = [1:3];
pl.time2baseline = [1];
pl.freqrange=[-0.1 0.1];
% pl.freqrange=[0 0];
pl.sub2plot = 1:numel(F.Subs2use);

pl.data_ind = []; pl.data_evo = [];
pl.RDKlabel = {'RDK1';'RDK2';'RDK3';'RDK4'};
pl.RDKposlabel2 = {'center';'center';'peri';'peri'};
pl.timelabel = cellfun(@(x) vararg2str(x),FFTdat.ffttimewin(pl.time2plot),'UniformOutput',false);

% extract data
% ([RDK num], con, time, sub)
pl.RDK.amp_data_ind = nan(numel(pl.RDKlabel),numel(F.conlabel_att),numel(pl.time2plot),numel(pl.sub2plot));
pl.RDK.amp_data_evo = pl.RDK.amp_data_ind;
pl.RDK.PLV_data_evo = pl.RDK.amp_data_ind;
pl.RDK.mscohere_data_evo = pl.RDK.amp_data_ind;



% the ones that are fixed
pl.RDK.con = permute(repmat((1:numel(F.conlabel_att))',[1 size(pl.RDK.amp_data_ind,[1 3 4])]),[2 1 3 4]);
pl.RDK.conlabel = permute(repmat((F.conlabel_att),[1 size(pl.RDK.amp_data_ind,[1 3 4])]),[2 1 3 4]);
pl.RDK.conflicker = permute(repmat((F.conlabel_flicker),[1 size(pl.RDK.amp_data_ind,[1 3 4])]),[2 1 3 4]);
pl.RDK.RDK_num = repmat(pl.RDKlabel,[1 size(pl.RDK.amp_data_ind, [2 3 4 ])]);
pl.RDK.RDK_id = repmat(F.con_RDKpresented',[1 1 size(pl.RDK.amp_data_ind, [3 4 ])]);
pl.RDK.RDK_isattended = permute(repmat(F.conRDKattended_label,[1 1 size(pl.RDK.amp_data_ind,[3 4])]),[2 1 3 4]);
pl.RDK.timewin = permute(repmat(pl.timelabel,[1 size(pl.RDK.amp_data_ind,[1 2 4])]), [2 3 1 4]);
pl.RDK.sub = permute(repmat(F.Subs2use(pl.sub2plot)',[1 size(pl.RDK.amp_data_ind,[1 2 3])]),[2 3 4 1]);
pl.RDK.RDK_pos2 = repmat(pl.RDKposlabel2,[1 size(pl.RDK.amp_data_ind, [2 3 4 ])]);

% the onses that change
pl.RDK.RDK_freq = pl.RDK.amp_data_ind;
pl.RDK.RDK_color = repmat({''},size(pl.RDK.amp_data_ind));
pl.RDK.RDK_electrodes = repmat({''},size(pl.RDK.amp_data_ind));

for i_sub = 1:numel(pl.sub2plot)
    for i_con = 1:numel(F.conlabel_att)
        for i_rdk = 1:numel(pl.RDKlabel)
            % index the actual RDK
            t.rdkidx = pl.RDK.RDK_id(i_rdk,i_con,1,1);
            if FFTdat.RDK(pl.sub2plot(i_sub)).RDK(t.rdkidx).in_flag == 1
                t.freq = FFTdat.RDK(pl.sub2plot(i_sub)).RDK(t.rdkidx).in_freq; % which SSVEP frequency?
            else
                t.freq = FFTdat.RDK(pl.sub2plot(i_sub)).RDK(t.rdkidx).freq; % which SSVEP frequency?
            end
            t.fidx = dsearchn(FFTdat.fftfreqs',(pl.freqrange+t.freq)');
            % index the electrode cluster
            t.elecidx = pl.elec2plot_i{pl.elec4conXrdk(i_con,i_rdk)};

            % extract data
            % induced
            pl.RDK.amp_data_ind(i_rdk,i_con,:,i_sub) = squeeze(mean( ...
                FFTdat.amp_data_ind(t.fidx(1):t.fidx(2),pl.elec2plot_i{1},i_con,pl.time2plot,pl.sub2plot(i_sub)),[1,2] ... %induced
                ));
            % evoked
            pl.RDK.amp_data_evo(i_rdk,i_con,:,i_sub) = squeeze(mean( ...
                FFTdat.amp_data_evo(t.fidx(1):t.fidx(2),pl.elec2plot_i{1},i_con,pl.time2plot,pl.sub2plot(i_sub)),[1,2] ... %induced
                ));
            % PLV
            pl.RDK.PLV_data_evo(i_rdk,i_con,:,i_sub) = squeeze(mean( ...
                FFTdat.PLV_data_evo(t.fidx(1):t.fidx(2),pl.elec2plot_i{1},i_con,pl.time2plot,pl.sub2plot(i_sub)),[1,2] ... %induced
                ));
            % mscohere
            pl.RDK.mscohere_data_evo(i_rdk,i_con,:,i_sub) = squeeze(mean( ...
                FFTdat.mscohere_data_evo(t.fidx(1):t.fidx(2),pl.elec2plot_i{1},i_con,pl.time2plot,pl.sub2plot(i_sub)),[1,2] ... %induced
                ));



            % write some data
            pl.RDK.RDK_freq(i_rdk,i_con,:,i_sub) = t.freq;
            pl.RDK.RDK_color(i_rdk,i_con,:,i_sub) = {FFTdat.RDK(pl.sub2plot(i_sub)).RDK(t.rdkidx).colnames};
            pl.RDK.RDK_electrodes(i_rdk,i_con,:,i_sub) = {vararg2str(pl.elec2plot{pl.elec4conXrdk(i_con,i_rdk)})};
        end
    end
end

% baseline corrected data
pl.RDK.amp_data_ind_bc = 100*(bsxfun(@rdivide, pl.RDK.amp_data_ind, pl.RDK.amp_data_ind(:,:,1,:))-1);
pl.RDK.amp_data_ind_bc_sub = bsxfun(@minus,  pl.RDK.amp_data_ind, pl.RDK.amp_data_ind(:,:,1,:));

pl.RDK.amp_data_evo_bc = 100*(bsxfun(@rdivide, pl.RDK.amp_data_evo, pl.RDK.amp_data_evo(:,:,1,:))-1);
pl.RDK.amp_data_evo_bc_sub = bsxfun(@minus, pl.RDK.amp_data_evo, pl.RDK.amp_data_evo(:,:,1,:));

pl.RDK.mscohere_data_evo_bc = 100*(bsxfun(@rdivide, pl.RDK.mscohere_data_evo, pl.RDK.mscohere_data_evo(:,:,1,:))-1);
pl.RDK.mscohere_data_evo_bc_sub = bsxfun(@minus, pl.RDK.mscohere_data_evo, pl.RDK.mscohere_data_evo(:,:,1,:));

pl.RDK.PLV_data_evo_bc = 100*(bsxfun(@rdivide, pl.RDK.PLV_data_evo, pl.RDK.PLV_data_evo(:,:,1,:))-1);
pl.RDK.PLV_data_evo_bc_sub = bsxfun(@minus, pl.RDK.PLV_data_evo, pl.RDK.PLV_data_evo(:,:,1,:));



% % trouble shooting do subtraction and modulation correspond?
% % pl.tdata = reshape(squeeze(pl.RDK.data_evo_bc(1:2,:,2,:)),[],size(pl.RDK.data_evo_bc,4));
% % pl.tdata(:,:,2) = reshape(squeeze(pl.RDK.data_evo_bc_sub(1:2,:,2,:)),[],size(pl.RDK.data_evo_bc,4)).*10;
% pl.tdata = reshape(squeeze(mean(pl.RDK.data_evo_bc(1:2,:,2,:),1)),[],size(pl.RDK.data_evo_bc,4));
% pl.tdata(:,:,2) = reshape(squeeze(mean(pl.RDK.data_evo_bc_sub(1:2,:,2,:),1)),[],size(pl.RDK.data_evo_bc,4)).*10;
% figure;
% for i_sub = 1:size(pl.tdata,2)
% %     plot([-0.25 +0.25]+i_sub, squeeze(pl.tdata(:,i_sub,:)),'Color',[0.3 0.3 0.3 0.5])
%     plot([-0.25 +0.25]+i_sub, sign(squeeze(pl.tdata(:,i_sub,:))),'Color',[0.3 0.3 0.3 0.5])
%     hold on
% end
% 
% pl.tdata = reshape(squeeze(diff(sign(pl.RDK.data_evo_bc(1:2,:,2,:)),1,1)),[],size(pl.RDK.data_evo_bc,4));
% pl.tdata(:,:,2) = reshape(squeeze(diff(sign(pl.RDK.data_evo_bc_sub(1:2,:,2,:)),1,1)),[],size(pl.RDK.data_evo_bc,4));
% figure;
% for i_sub = 1:size(pl.tdata,2)
%     plot([-0.25 +0.25]+i_sub, ...
%         squeeze(pl.tdata(:,i_sub,:))+repmat(((randn(size(pl.tdata,1),1)-0.5).*0.1),1,2), ...
%         'Color',[0.3 0.3 0.3 0.5])
% %     plot([-0.25 +0.25]+i_sub, sign(squeeze(pl.tdata(:,i_sub,:))),'Color',[0.3 0.3 0.3 0.5])
%     hold on
% end
% figure; plot(sign(pl.RDK.data_evo_bc_sub(:)), sign(pl.RDK.data_evo_bc(:)))
% figure; plot(sign(R_Mat.all_table.modulation_evoked), sign(R_Mat.all_table.subtraction_evoked))


% % do some interim plotting for checking everything
% t.tdata = cat(3, ...
%     [squeeze(mean(pl.RDK.data_evo(1,[1 3 5],2,:),2)) squeeze(mean(pl.RDK.data_evo(1,[2 4 6],2,:),2))], ...
%     [squeeze(mean(pl.RDK.data_evo(2,[2 4 6],2,:),2)) squeeze(mean(pl.RDK.data_evo(2,[1 3 5],2,:),2))]);
% figure; boxplot(mean(t.tdata,3))
% t.tidx = strcmp(pl.RDK.RDK_pos1,'center') & strcmp(pl.RDK.timewin,'[0.5 1.5] ') & pl.RDK.RDK_freq == 29 ;
% t.tdata_ev0 =  pl.RDK.data_evo_bc;  t.tdata_ev0(~t.tidx)= nan;
% t.tdata = cat(3, ...
%     [squeeze(mean(t.tdata_ev0(1,[1 3 5],2,:),2,'omitnan')) squeeze(mean(t.tdata_ev0(1,[2 4 6],2,:),2,'omitnan'))], ...
%     [squeeze(mean(t.tdata_ev0(2,[2 4 6],2,:),2,'omitnan')) squeeze(mean(t.tdata_ev0(2,[1 3 5],2,:),2,'omitnan'))]);
% figure; boxplot(mean(t.tdata,3,'omitnan'))
% figure; plot(mean(t.tdata,3,'omitnan')')
% t.tdata2 =diff(mean(t.tdata,3,'omitnan'),1,2);
% figure; boxplot(diff(mean(t.tdata,3,'omitnan'),1,2));
% figure; histfit(diff(mean(t.tdata,3,'omitnan'),1,2),20)
% [tt.h,tt.p,tt.ci,tt.stats] = ttest(diff(mean(t.tdata,3,'omitnan'),1,2));

R_Mat.all = [{'amplitude_induced','amplitude_evoked','mscohere','PLV', ...
    'modulation_amp_induced','modulation_amp_evoked','modulation_mscohere','modulation_PLV', ...
    'subtraction_amp_induced','subtraction_amp_evoked','subtraction_mscohere','subtraction_PLV', ...
    'subjects', 'condition', 'con_flicker','time', ...
    'RDK_id', 'RDK_num', 'RDK_position', 'RDK_freq', 'RDK_color', 'RDK_isattended', 'RDK_electrodes'}; ...
    num2cell([pl.RDK.amp_data_ind(:) pl.RDK.amp_data_evo(:) pl.RDK.mscohere_data_evo(:) pl.RDK.PLV_data_evo(:) ...
    pl.RDK.amp_data_ind_bc(:) pl.RDK.amp_data_evo_bc(:) pl.RDK.mscohere_data_evo_bc(:) pl.RDK.PLV_data_evo_bc(:) ...
    pl.RDK.amp_data_ind_bc_sub(:) pl.RDK.amp_data_evo_bc_sub(:) pl.RDK.mscohere_data_evo_bc_sub(:) pl.RDK.PLV_data_evo_bc_sub(:) ...
    pl.RDK.sub(:) pl.RDK.con(:)]) pl.RDK.conflicker(:) pl.RDK.timewin(:) ...
    num2cell(pl.RDK.RDK_id(:)) pl.RDK.RDK_num(:) pl.RDK.RDK_pos2(:) num2cell(pl.RDK.RDK_freq(:)) pl.RDK.RDK_color(:) ...
    pl.RDK.RDK_isattended(:) pl.RDK.RDK_electrodes(:)
    ];

R_Mat.all_table = cell2table(R_Mat.all(2:end,:), "VariableNames",R_Mat.all(1,:));


t.path = 'C:\Dropboxdata\Dropbox\work\R-statistics\experiments\ssvep_fshifthighfreqflick\data_in';
% t.path = 'C:\Users\EEG\Documents\R\Christopher\analysis_R_ssvep_fshift_perirr\data_in';
t.datestr = datestr(now,'mm-dd-yyyy_HH-MM');
% write to textfile

% writetable(R_Mat.all_table,fullfile(t.path,sprintf('FFT_Amp_data_largeclusters_%s.csv',t.datestr)),'Delimiter',';')


