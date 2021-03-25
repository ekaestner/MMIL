clear; clc;
addpath /home/ekaestne/fieldtrip-20131031/;
ft_defaults

semantic = 1;
repeated = 0;

lfp = 1;
hgp = 0;

subjn = 1;

% chan = {'2>G15' '2>G21' '2>G22' '2>G23' '2>G27' '2>G28' '2>G29' '2>G30' '2>G31'  '2>G49' '2>G57' ...
%     '2>G58' '2>G63' '2>G64' '2>SF5' '2>SF5' '2>PC8' '2>MT3' '2>FT3' '2>FT4' '2>PT3' '2>TO1' '2>TO2'}; %NY007 HGP repeated
% chan = {'2>G22' '2>G28'};  %NY007 HGP semantic
% chan = {'2>G8' '2>G14' '2>G15' '2>G16' '2>G22' '2>G27' '2>G31' '2>G57' '2>G51' '2>AT2' '2>MT3' '2>FT3' ...
%     '2>FT4' '2>D3' '2>D4' '2>D5' '2>D6' '2>D7' '2>D2(4)' '2>D2(5)'}; % NY007 LFP repeated
chan = {'2>G22' '2>G29' '2>G30' '2>G58' '2>MT3' '2>PT3'}; % NY007 LFP semantic

% chan = {'2>G3' '2>G4' '2>G7' '2>G12' '2>G21' '2>G28' '2>G36' '2>G37' '2>G37'  '2>G39' '2>G42' ...
%      '2>G43' '2>G44' '2>G47' '2>G48' '2>G56' '2>AT3' '2>MT6' '2>PT1' '2>TO8' '2>EPT1'}; %NY008 HGP repeated
% chan = {'2>G4' '2>G6' '2>G7' '2>G8' '2>G10' '2>G13' '2>G19' '2>G21' '2>G22' '2>G24' '2>G29' '2>G39' ... 
%     '2>AT1' '2>MT2' '2>PT2' '2>PT4' '2>EPT4' '2>EPT5'}; % NY008 repeated
% chan = {'2>EPT1'}; % NY008 Semantic

% chan = {'TO3' 'TO4' 'TO5' 'TO6' 'TO7' 'PT3' 'PT4' 'PT6'}; %NY011 HGP repeated
% chan = {'G22' 'G24' 'G28' 'G29' 'G31' 'G34' 'G47' 'G48' 'G50' 'G56' 'G63' 'AT1' 'AT2' ...
%     'MT2' 'MT3' 'MT4' 'MT6' 'DB1'}; %NY011 LFP repeated
% chan = {'PT6' 'DC1' 'DC2' 'DC3' 'DC4' 'DD1' 'DD2' 'DD3' 'DD4' 'DD5'}; %NY011 LFP semantic

% chan = {'G6' 'G7' 'G36' 'G37' 'G38' 'G46' 'G55' 'G56' 'PT1' 'PT2' 'P1'}; %NY013 HGP repeated
% chan = {'G7' 'G56' 'PT1' 'AT1' 'AT2' 'PG13'}; %NY013 LFP repeated
% chan = {'G7'}; %NY013 LFP semantic

% chan = {'2>G27' '2>G30' '2>G40' '2>G55' '2>G51' '2>G59'}; %NY014 HGP repeated
% chan = {'2>G41' '2>G61'}; %NY014 LFP repeated

% chan = {'G5' 'G8' 'G13' 'G15' 'G21' 'G22' 'G23' 'G52' 'G55' 'PMT5' 'PT1'}; %NY018 HGP repeated
% chan = {'G14'}; %NY018 LFP repeated
% chan = {'G21'}; %NY018 LFP semantic

% chan = {'G6' 'G14' 'G19' 'G20' 'G16' 'G21' 'G22' 'G23' 'G26' 'G27' 'G28' 'G29' 'G30' 'G31' 'PT3' 'PT1' ...
%     'O2' 'O4' 'O5' 'PF5'}; %NY023 HGP repeated
% chan = {'G21' 'G27' 'G29'}; %NY023 HGP semantic
% chan = {'PT5' 'PF5'}; %NY023 LFP repeated
% chan = {'G20' 'G8' 'G28' 'O4' 'O5'};

% chan = {'G12' 'G13' 'G21' 'G32' 'G39' 'G43' 'G56' 'G59' 'G60' 'G61' 'G64' 'PT3' 'PT4'}; %NY024 HGP repeated
% chan = {'G5' 'G7' 'G13' 'G16' 'G32' 'G36' 'G38' 'G39' 'G44' 'G64' 'F3' 'F4' 'AP1' 'AP2' 'AT2' ...
%     'DSA7' 'DSA8' 'DSA1' 'DSA2' 'DSA3'}; %NY024 LFP repeated
% chan = {'G8' 'G16' 'G32' '05' 'DSA6'}; %NY024 LFP semantic


subj = {'NY007' ...
    'NY008' ...
    'NY011' ...
    'NY013' ...
    'NY014' ...
    'NY018' ...
    'NY023' ...
    'NY024' ...
    };
subj = subj{subjn};

indir = '/space/mdeh3/9/halgdev/projects/mmilanguage/iSASZ/proc_data';
fls = subdir(fullfile(indir,'*proc.mat'));
flsn = find(~cellfun(@isempty,strfind({fls.name},subj)));

load(fls(flsn(1)).name)
load(fls(flsn(2)).name)

if lfp
    nchan_lfp = cell2mat(cellfun(@find,cellfun(@(x) strcmpi(epoch_data_lfp_preproc.label,x),chan,'UniformOutput',0),'UniformOutput',0));
end
if hgp
    nchan_hgp = cell2mat(cellfun(@find,cellfun(@(x) strcmpi(epoch_data_hgp_preproc.label,x),chan,'UniformOutput',0),'UniformOutput',0));
end

% Preprocess Data
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [0.5 15];
epoch_data_lfp_preproc = ft_func(@ft_preprocessing,cfg,epoch_data_lfp_preproc);
epoch_data_hgp_preproc = ft_func(@ft_preprocessing,cfg,epoch_data_hgp_preproc);

cfg = [];
cfg.toilim = [-0.3 1.2];
cfg.trials = 'all';
epoch_data_lfp_preproc = ft_func(@ft_redefinetrial,cfg,epoch_data_lfp_preproc);
epoch_data_hgp_preproc = ft_func(@ft_redefinetrial,cfg,epoch_data_hgp_preproc);

if ~exist([indir '/' subj '/' 'vis'],'dir')
    mkdir([indir '/' subj '/' 'vis'])
end

for ic = 1:length(chan)
    
    %% Repeated
    if repeated && lfp
        
        load([indir '/' subj '/' 'stat_data_sa_repeated_lfp.mat']) % stat_data_sa_repeated
        load([indir '/' subj '/' 'stat_data_sz_repeated_lfp.mat']) % stat_data_sz_repeated
        
        cfg = [];
        cfg.return_events = 0;
        cfg.old_events = {[1 2 3 4] [5 6 7 8] [11 12 13 14] [15 16 17 18]};
        cfg.new_events = {1 2 3 4};
        epoch_data_lfp_preproc = ft_redefine_events(cfg,epoch_data_lfp_preproc);
        
        ft_dat_lfp = cat(3,epoch_data_lfp_preproc.trial{:});
        
        evs = [1 2 3 4];
        for iev = 1:length(evs)
            ev_ind_lfp{iev} = find(epoch_data_lfp_preproc.trialinfo==evs(iev));
        end
        
        maxy = zeros(1,length(ev_ind_lfp));
        miny = zeros(1,length(ev_ind_lfp));
        
%         col = {'c' 'y' 'c' 'y'};
%         col_dev = {[0.3 0.7 0.7] [0.7 0.7 0.3] [0.3 0.7 0.7] [0.7 0.7 0.3]};
        col = {[0.3 0.7 0.7] [0.7 0.7 0.3] [0.3 0.7 0.7] [0.7 0.7 0.3]};
        
        %% LFP       
        figure('Visible','Off');
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 5])
        
        ut_hand(1) = subplot(1,2,1);
        for iev = 3:4
            hold on
            
            mdat = squeeze(ft_dat_lfp(nchan_lfp(ic),:,ev_ind_lfp{iev}))';
            ev_mean = mean(mdat);
            
            plot(epoch_data_lfp_preproc.time{1},ev_mean,'Color',col{iev},'LineWidth',2);
            
%             % Deviation
            lowb  = ev_mean-(std(mdat)/sqrt(size(mdat,1)));
            highb = ev_mean+(std(mdat)/sqrt(size(mdat,1)));
%             
%             % Add measures of variance
%             patch([epoch_data_lfp_preproc.time{1} epoch_data_lfp_preproc.time{1}(end:-1:1)],[highb lowb(end:-1:1)],col_dev{iev})
%             chldrn = get(gca,'Children');
%             lnchldrn = findobj(gca,'Type','line');
%             [~,iH] = sort(ismember(chldrn,lnchldrn),1,'descend');
%             set(gca,'Children',chldrn(iH));
            
            % For plotting purposes to show the highest value
            maxy(1,iev) = max(highb);
            miny(1,iev) = min(lowb);
        end
        
        % Stats
        maxy_curr = max(maxy(:));
        miny_curr = min(miny(:));
        
        stats = stat_data_sa_repeated.stats(1).prob(nchan_lfp(ic),:) < 0.05;
        
        if stats(1) == 1; stats(1) = 0; end
        if stats(end) == 1; stats(end) = 0; end
        sig_periods = ts_crossing(stats(:),[],0.5); % - EJK change to 'crossing.m'
        sig_1st = sig_periods(1:2:end);
        sig_2nd = sig_periods(2:2:end);
        
        lnchldrn = findobj(gca,'Type','line');
        patchchldrn = findobj(gca,'Type','patch');
        combchldrn = [lnchldrn ; patchchldrn];
        
        for isig = 1:length(sig_1st)
            ftime = epoch_data_lfp_preproc.time{1}(sig_1st(isig):sig_2nd(isig)) ; btime = epoch_data_lfp_preproc.time{1}(sig_2nd(isig):-1:sig_1st(isig));
            patch([ftime btime],[repmat(maxy_curr,1,sig_2nd(isig)-sig_1st(isig)+1) repmat(miny_curr,1,sig_2nd(isig)-sig_1st(isig)+1)],[0.7 0.7 0.7])
            chldrn = get(gca,'Children');
            [~,iH] = sort(ismember(chldrn,combchldrn),1,'descend');
            set(gca,'Children',chldrn(iH));
        end
        
        hold off
        
        ut_hand(2) = subplot(1,2,2);
        for iev = 1:2
            
            hold on
            
            mdat = squeeze(ft_dat_lfp(nchan_lfp(ic),:,ev_ind_lfp{iev}))';
            ev_mean = mean(mdat);
            
            plot(epoch_data_lfp_preproc.time{1},ev_mean,'Color',col{iev},'LineWidth',2);
            
%             % Deviation
            lowb  = ev_mean-(std(mdat)/sqrt(size(mdat,1)));
            highb = ev_mean+(std(mdat)/sqrt(size(mdat,1)));
%             
%             % Add measures of variance
%             patch([epoch_data_lfp_preproc.time{1} epoch_data_lfp_preproc.time{1}(end:-1:1)],[highb lowb(end:-1:1)],col_dev{iev})
%             chldrn = get(gca,'Children');
%             lnchldrn = findobj(gca,'Type','line');
%             [~,iH] = sort(ismember(chldrn,lnchldrn),1,'descend');
%             set(gca,'Children',chldrn(iH));
            
            % For plotting purposes to show the highest value
            maxy(1,iev) = max(highb);
            miny(1,iev) = min(lowb);
        end
        
        % Stats
        maxy_curr = max(maxy(:));
        miny_curr = min(miny(:));
        
        stats = stat_data_sz_repeated.stats(1).prob(nchan_lfp(ic),:) < 0.05;
        
        if stats(1) == 1; stats(1) = 0; end
        if stats(end) == 1; stats(end) = 0; end
        sig_periods = ts_crossing(stats(:),[],0.5); % - EJK change to 'crossing.m'
        sig_1st = sig_periods(1:2:end);
        sig_2nd = sig_periods(2:2:end);
        
        lnchldrn = findobj(gca,'Type','line');
        patchchldrn = findobj(gca,'Type','patch');
        combchldrn = [lnchldrn ; patchchldrn];
        
        for isig = 1:length(sig_1st)
            ftime = epoch_data_lfp_preproc.time{1}(sig_1st(isig):sig_2nd(isig)) ; btime = epoch_data_lfp_preproc.time{1}(sig_2nd(isig):-1:sig_1st(isig));
            patch([ftime btime],[repmat(maxy_curr,1,sig_2nd(isig)-sig_1st(isig)+1) repmat(miny_curr,1,sig_2nd(isig)-sig_1st(isig)+1)],[0.7 0.7 0.7])
            chldrn = get(gca,'Children');
            [~,iH] = sort(ismember(chldrn,combchldrn),1,'descend');
            set(gca,'Children',chldrn(iH));
        end
        
        % Remove Axis Information
        set(gca,'YTickLabel',[])
        set(gca,'XTickLabel',[])
        
        linkaxes(ut_hand,'xy');
        xlim([epoch_data_lfp_preproc.time{1}(1) epoch_data_lfp_preproc.time{1}(end)])
        
        maxy = max(max(maxy));
        miny = min(min(miny));
        ylim([miny-miny*0.1 maxy+maxy*0.1])
        
        tick_y = (miny-miny*0.1):((abs(miny-miny*0.1)+abs(maxy+maxy*0.1))/3):(maxy+maxy*0.1);
        tick_x = epoch_data_lfp_preproc.time{1}(1):0.2:epoch_data_lfp_preproc.time{1}(end);
        
        set(ut_hand(:),'YTick',tick_y)
        set(ut_hand(:),'XTick',tick_x)
        
        linkaxes(ut_hand,'off');
        
        % Add stimulus line
        for il = 1:length(ut_hand)
            subplot(1,2,il);
            line('XData',[0 0],'YData',[floor(miny-miny*0.1) ceil(maxy+maxy*0.1)],'LineWidth',2,'Color','k');
        end
        
        labe = strrep(epoch_data_lfp_preproc.label{nchan_lfp(ic)},'2>',''); 
        print(gcf,[indir '/' subj '/' 'vis' '/' subj '_' labe '_lfp_repeated'],'-dpng','-r500')
        close all
        
        cfg = [];
        cfg.return_events = 1;
        epoch_data_lfp_preproc = ft_redefine_events(cfg,epoch_data_lfp_preproc);
                
    end
    
    %% HGP
    if repeated && hgp
        
        load([indir '/' subj '/' 'stat_data_sa_repeated_hgp.mat']) % stat_data_sa_repeated
        load([indir '/' subj '/' 'stat_data_sz_repeated_hgp.mat']) % stat_data_sz_repeated
        
        cfg = [];
        cfg.return_events = 0;
        cfg.old_events = {[1 2 3 4] [5 6 7 8] [11 12 13 14] [15 16 17 18]};
        cfg.new_events = {1 2 3 4};
        epoch_data_hgp_preproc = ft_redefine_events(cfg,epoch_data_hgp_preproc);
        
        ft_dat_hgp = cat(3,epoch_data_hgp_preproc.trial{:});
        
        load([indir '/' subj '/' 'stat_data_sa_repeated_hgp.mat']) % stat_data_sa_repeated
        load([indir '/' subj '/' 'stat_data_sz_repeated_hgp.mat']) % stat_data_sz_repeated
        
        evs = [1 2 3 4];
        for iev = 1:length(evs)
            ev_ind_hgp{iev} = find(epoch_data_hgp_preproc.trialinfo==evs(iev));
        end
        
        maxy = zeros(1,length(ev_ind_hgp));
        miny = zeros(1,length(ev_ind_hgp));
        
%         col = {'c' 'y' 'c' 'y'};
%         col_dev = {[0.3 0.7 0.7] [0.7 0.7 0.3] [0.3 0.7 0.7] [0.7 0.7 0.3]};
        col = {[0.3 0.7 0.7] [0.7 0.7 0.3] [0.3 0.7 0.7] [0.7 0.7 0.3]};
        
        % hgp
        %     figure()
        figure('Visible','Off');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 5])
        
        hold on
        
        ut_hand(1) = subplot(1,2,1);
        for iev = 3:4
            
            hold on
            
            mdat = squeeze(ft_dat_hgp(nchan_hgp(ic),:,ev_ind_hgp{iev}))';
            ev_mean = mean(mdat);
            
            plot(epoch_data_hgp_preproc.time{1},ev_mean,'Color',col{iev},'LineWidth',2);
            
%             % Deviation
            lowb  = ev_mean-(std(mdat)/sqrt(size(mdat,1)));
            highb = ev_mean+(std(mdat)/sqrt(size(mdat,1)));
%             
%             % Add measures of variance
%             patch([epoch_data_hgp_preproc.time{1} epoch_data_hgp_preproc.time{1}(end:-1:1)],[highb lowb(end:-1:1)],col_dev{iev})
%             chldrn = get(gca,'Children');
%             lnchldrn = findobj(gca,'Type','line');
%             [~,iH] = sort(ismember(chldrn,lnchldrn),1,'descend');
%             set(gca,'Children',chldrn(iH));
            
            % For plotting purposes to show the highest value
            maxy(1,iev) = max(highb);
            miny(1,iev) = min(lowb);
        end
        
        % Stats
        maxy_curr = max(maxy(:));
        miny_curr = min(miny(:));
        
        stats = stat_data_sa_repeated.stats(1).prob(nchan_hgp(ic),:) < 0.05;
        
        if stats(1) == 1; stats(1) = 0; end
        if stats(end) == 1; stats(end) = 0; end
        sig_periods = ts_crossing(stats(:),[],0.5); % - EJK change to 'crossing.m'
        sig_1st = sig_periods(1:2:end);
        sig_2nd = sig_periods(2:2:end);
        
        lnchldrn = findobj(gca,'Type','line');
        patchchldrn = findobj(gca,'Type','patch');
        combchldrn = [lnchldrn ; patchchldrn];
        
        for isig = 1:length(sig_1st)
            ftime = epoch_data_hgp_preproc.time{1}(sig_1st(isig):sig_2nd(isig)) ; btime = epoch_data_hgp_preproc.time{1}(sig_2nd(isig):-1:sig_1st(isig));
            patch([ftime btime],[repmat(maxy_curr,1,sig_2nd(isig)-sig_1st(isig)+1) repmat(miny_curr,1,sig_2nd(isig)-sig_1st(isig)+1)],[0.7 0.7 0.7])
            chldrn = get(gca,'Children');
            [~,iH] = sort(ismember(chldrn,combchldrn),1,'descend');
            set(gca,'Children',chldrn(iH));
        end
        
        % Remove Axis Information
        set(gca,'YTickLabel',[])
        set(gca,'XTickLabel',[])
        
        hold off
        
        ut_hand(2) = subplot(1,2,2);
        for iev = 1:2
            
            hold on
            
            mdat = squeeze(ft_dat_hgp(nchan_hgp(ic),:,ev_ind_hgp{iev}))';
            ev_mean = mean(mdat);
            
            plot(epoch_data_hgp_preproc.time{1},ev_mean,'Color',col{iev},'LineWidth',2);
            
%             % Deviation
            lowb  = ev_mean-(std(mdat)/sqrt(size(mdat,1)));
            highb = ev_mean+(std(mdat)/sqrt(size(mdat,1)));
%             
%             % Add measures of variance
%             patch([epoch_data_hgp_preproc.time{1} epoch_data_hgp_preproc.time{1}(end:-1:1)],[highb lowb(end:-1:1)],col_dev{iev})
%             chldrn = get(gca,'Children');
%             lnchldrn = findobj(gca,'Type','line');
%             [~,iH] = sort(ismember(chldrn,lnchldrn),1,'descend');
%             set(gca,'Children',chldrn(iH));
            
            % For plotting purposes to show the highest value
            maxy(1,iev) = max(highb);
            miny(1,iev) = min(lowb);
        end
        
        % Stats
        maxy_curr = max(maxy(:));
        miny_curr = min(miny(:));
        
        stats = stat_data_sz_repeated.stats(1).prob(nchan_hgp(ic),:) < 0.05;
        
        if stats(1) == 1; stats(1) = 0; end
        if stats(end) == 1; stats(end) = 0; end
        sig_periods = ts_crossing(stats(:),[],0.5); % - EJK change to 'crossing.m'
        sig_1st = sig_periods(1:2:end);
        sig_2nd = sig_periods(2:2:end);
        
        lnchldrn = findobj(gca,'Type','line');
        patchchldrn = findobj(gca,'Type','patch');
        combchldrn = [lnchldrn ; patchchldrn];
        
        for isig = 1:length(sig_1st)
            ftime = epoch_data_hgp_preproc.time{1}(sig_1st(isig):sig_2nd(isig)) ; btime = epoch_data_hgp_preproc.time{1}(sig_2nd(isig):-1:sig_1st(isig));
            patch([ftime btime],[repmat(maxy_curr,1,sig_2nd(isig)-sig_1st(isig)+1) repmat(miny_curr,1,sig_2nd(isig)-sig_1st(isig)+1)],[0.7 0.7 0.7])
            chldrn = get(gca,'Children');
            [~,iH] = sort(ismember(chldrn,combchldrn),1,'descend');
            set(gca,'Children',chldrn(iH));
        end
        
        % Remove Axis Information
        set(gca,'YTickLabel',[])
        set(gca,'XTickLabel',[])
        
        linkaxes(ut_hand,'xy');
        xlim([epoch_data_hgp_preproc.time{1}(1) epoch_data_hgp_preproc.time{1}(end)])
        
        maxy = max(max(maxy));
        miny = min(min(miny));
        ylim([miny-miny*0.1 maxy+maxy*0.1])
        
        tick_y = (miny-miny*0.1):((abs(miny-miny*0.1)+abs(maxy+maxy*0.1))/3):(maxy+maxy*0.1);
        tick_x = epoch_data_hgp_preproc.time{1}(1):0.2:epoch_data_hgp_preproc.time{1}(end);
        
        set(ut_hand(:),'YTick',tick_y)
        set(ut_hand(:),'XTick',tick_x)
        
        linkaxes(ut_hand,'off');
        
        % Label first left axis
        set(ut_hand(1),'YTickLabel',tick_y,'FontSize',21)
        set(ut_hand(1),'XTickLabel',tick_x,'FontSize',21)
        
        % Add stimulus line
        for il = 1:length(ut_hand)
            subplot(1,2,il);
            line('XData',[0 0],'YData',[floor(miny-miny*0.1) ceil(maxy+maxy*0.1)],'LineWidth',2,'Color','k');
        end
        
        labe = strrep(epoch_data_lfp_preproc.label{nchan_lfp(ic)},'2>',''); 
        print(gcf,[indir '/' subj '/' 'vis' '/' subj '_' labe '_hgp_repeated'],'-dpng','-r500')
        close all
        
        cfg = [];
        cfg.return_events = 1;
        epoch_data_hgp_preproc = ft_redefine_events(cfg,epoch_data_hgp_preproc);
        
    end
    
end

%% Semantic
for ic = 1:length(nchan_lfp)
    
    if semantic && lfp
        %% semantic
        load([indir '/' subj '/' 'stat_data_sa_semantic_lfp.mat']) % stat_data_sa_semantic
        load([indir '/' subj '/' 'stat_data_sz_semantic_lfp.mat']) % stat_data_sa_semantic
        
        cfg = [];
        cfg.return_events = 0;
        cfg.old_events = {[1 3 5 7] [2 4 6 8] [11 13 15 17] [12 14 16 18]};
        cfg.new_events = {1 2 3 4};
        epoch_data_lfp_preproc = ft_redefine_events(cfg,epoch_data_lfp_preproc);
        
        ft_dat_lfp = cat(3,epoch_data_lfp_preproc.trial{:});
        
        evs = [1 2 3 4];
        for iev = 1:length(evs)
            ev_ind_lfp{iev} = find(epoch_data_lfp_preproc.trialinfo==evs(iev));
        end
        
        maxy = zeros(1,length(ev_ind_lfp));
        miny = zeros(1,length(ev_ind_lfp));
        
%         col = {'g' 'm' 'g' 'm'};
%         col_dev = {[0.3 0.7 0.3] [0.7 0.3 0.7] [0.3 0.7 0.3] [0.7 0.3 0.7]};
        col = {[0.3 0.7 0.3] [0.7 0.3 0.7] [0.3 0.7 0.3] [0.7 0.3 0.7]};
        
        %% LFP
        load([indir '/' subj '/' 'stat_data_sa_semantic_lfp.mat']) % stat_data_sa_semantic
        load([indir '/' subj '/' 'stat_data_sz_semantic_lfp.mat']) % stat_data_sz_semantic
        
        figure('Visible','Off');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 5])
        
        ut_hand(1) = subplot(1,2,1);
        for iev = 3:4
            hold on
            
            mdat = squeeze(ft_dat_lfp(nchan_lfp(ic),:,ev_ind_lfp{iev}))';
            ev_mean = mean(mdat);
            
            plot(epoch_data_lfp_preproc.time{1},ev_mean,'Color',col{iev},'LineWidth',2);
            
%             % Deviation
            lowb  = ev_mean-(std(mdat)/sqrt(size(mdat,1)));
            highb = ev_mean+(std(mdat)/sqrt(size(mdat,1)));
%             
%             % Add measures of variance
%             patch([epoch_data_lfp_preproc.time{1} epoch_data_lfp_preproc.time{1}(end:-1:1)],[highb lowb(end:-1:1)],col_dev{iev})
%             chldrn = get(gca,'Children');
%             lnchldrn = findobj(gca,'Type','line');
%             [~,iH] = sort(ismember(chldrn,lnchldrn),1,'descend');
%             set(gca,'Children',chldrn(iH));
            
            % For plotting purposes to show the highest value
            maxy(1,iev) = max(highb);
            miny(1,iev) = min(lowb);
        end
        
        % Stats
        maxy_curr = max(maxy(:));
        miny_curr = min(miny(:));
        
        stats = stat_data_sa_semantic.stats(1).prob(nchan_lfp(ic),:) < 0.05;
        
        if stats(1) == 1; stats(1) = 0; end
        if stats(end) == 1; stats(end) = 0; end
        sig_periods = ts_crossing(stats(:),[],0.5); % - EJK change to 'crossing.m'
        sig_1st = sig_periods(1:2:end);
        sig_2nd = sig_periods(2:2:end);
        
        lnchldrn = findobj(gca,'Type','line');
        patchchldrn = findobj(gca,'Type','patch');
        combchldrn = [lnchldrn ; patchchldrn];
        
        for isig = 1:length(sig_1st)
            ftime = epoch_data_lfp_preproc.time{1}(sig_1st(isig):sig_2nd(isig)) ; btime = epoch_data_lfp_preproc.time{1}(sig_2nd(isig):-1:sig_1st(isig));
            patch([ftime btime],[repmat(maxy_curr,1,sig_2nd(isig)-sig_1st(isig)+1) repmat(miny_curr,1,sig_2nd(isig)-sig_1st(isig)+1)],[0.7 0.7 0.7])
            chldrn = get(gca,'Children');
            [~,iH] = sort(ismember(chldrn,combchldrn),1,'descend');
            set(gca,'Children',chldrn(iH));
        end
        
        hold off
        
        ut_hand(2) = subplot(1,2,2);
        for iev = 1:2
            
            hold on
            
            mdat = squeeze(ft_dat_lfp(nchan_lfp(ic),:,ev_ind_lfp{iev}))';
            ev_mean = mean(mdat);
            
            plot(epoch_data_lfp_preproc.time{1},ev_mean,'Color',col{iev},'LineWidth',2);
            
%             % Deviation
            lowb  = ev_mean-(std(mdat)/sqrt(size(mdat,1)));
            highb = ev_mean+(std(mdat)/sqrt(size(mdat,1)));
%             
%             % Add measures of variance
%             patch([epoch_data_lfp_preproc.time{1} epoch_data_lfp_preproc.time{1}(end:-1:1)],[highb lowb(end:-1:1)],col_dev{iev})
%             chldrn = get(gca,'Children');
%             lnchldrn = findobj(gca,'Type','line');
%             [~,iH] = sort(ismember(chldrn,lnchldrn),1,'descend');
%             set(gca,'Children',chldrn(iH));
            
            % For plotting purposes to show the highest value
            maxy(1,iev) = max(highb);
            miny(1,iev) = min(lowb);
        end
        
        % Stats
        maxy_curr = max(maxy(:));
        miny_curr = min(miny(:));
        
        stats = stat_data_sz_semantic.stats(1).prob(nchan_lfp(ic),:) < 0.05;
        
        if stats(1) == 1; stats(1) = 0; end
        if stats(end) == 1; stats(end) = 0; end
        sig_periods = ts_crossing(stats(:),[],0.5); % - EJK change to 'crossing.m'
        sig_1st = sig_periods(1:2:end);
        sig_2nd = sig_periods(2:2:end);
        
        lnchldrn = findobj(gca,'Type','line');
        patchchldrn = findobj(gca,'Type','patch');
        combchldrn = [lnchldrn ; patchchldrn];
        
        for isig = 1:length(sig_1st)
            ftime = epoch_data_lfp_preproc.time{1}(sig_1st(isig):sig_2nd(isig)) ; btime = epoch_data_lfp_preproc.time{1}(sig_2nd(isig):-1:sig_1st(isig));
            patch([ftime btime],[repmat(maxy_curr,1,sig_2nd(isig)-sig_1st(isig)+1) repmat(miny_curr,1,sig_2nd(isig)-sig_1st(isig)+1)],[0.7 0.7 0.7])
            chldrn = get(gca,'Children');
            [~,iH] = sort(ismember(chldrn,combchldrn),1,'descend');
            set(gca,'Children',chldrn(iH));
        end
        
        % Remove Axis Information
        set(gca,'YTickLabel',[])
        set(gca,'XTickLabel',[])
        
        linkaxes(ut_hand,'xy');
        xlim([epoch_data_lfp_preproc.time{1}(1) epoch_data_lfp_preproc.time{1}(end)])
        
        maxy = max(max(maxy));
        miny = min(min(miny));
        ylim([miny-miny*0.1 maxy+maxy*0.1])
        
        tick_y = (miny-miny*0.1):((abs(miny-miny*0.1)+abs(maxy+maxy*0.1))/3):(maxy+maxy*0.1);
        tick_x = epoch_data_lfp_preproc.time{1}(1):0.2:epoch_data_lfp_preproc.time{1}(end);
        
        set(ut_hand(:),'YTick',tick_y)
        set(ut_hand(:),'XTick',tick_x)
        
        linkaxes(ut_hand,'off');
        
        % Add stimulus line
        for il = 1:length(ut_hand)
            subplot(1,2,il);
            line('XData',[0 0],'YData',[floor(miny-miny*0.1) ceil(maxy+maxy*0.1)],'LineWidth',2,'Color','k');
        end
        
         labe = strrep(epoch_data_lfp_preproc.label{nchan_lfp(ic)},'2>',''); 
        print(gcf,[indir '/' subj '/' 'vis' '/' subj '_' labe '_lfp_semantic'],'-dpng','-r500')
        close all
        
        cfg = [];
        cfg.return_events = 1;
        epoch_data_lfp_preproc = ft_redefine_events(cfg,epoch_data_lfp_preproc);
   
    end
    
    %% HGP
    if semantic && hgp
        
        load([indir '/' subj '/' 'stat_data_sa_semantic_hgp.mat']) % stat_data_sa_semantic
        load([indir '/' subj '/' 'stat_data_sz_semantic_hgp.mat']) % stat_data_sa_semantic
        
        cfg = [];
        cfg.return_events = 0;
        cfg.old_events = {[1 3 5 7] [2 4 6 8] [11 13 15 17] [12 14 16 18]};
        cfg.new_events = {1 2 3 4};
        epoch_data_hgp_preproc = ft_redefine_events(cfg,epoch_data_hgp_preproc);
        
        ft_dat_hgp = cat(3,epoch_data_hgp_preproc.trial{:});
        
        load([indir '/' subj '/' 'stat_data_sa_semantic_hgp.mat']) % stat_data_sa_semantic
        load([indir '/' subj '/' 'stat_data_sz_semantic_hgp.mat']) % stat_data_sz_semantic
        
        evs = [1 2 3 4];
        for iev = 1:length(evs)
            ev_ind_hgp{iev} = find(epoch_data_hgp_preproc.trialinfo==evs(iev));
        end
        
        maxy = zeros(1,length(ev_ind_hgp));
        miny = zeros(1,length(ev_ind_hgp));
        
%         col = {'g' 'm' 'g' 'm'};
%         col_dev = {[0.3 0.7 0.3] [0.7 0.3 0.7] [0.3 0.7 0.3] [0.7 0.3 0.7]};
        col = {[0.3 0.7 0.3] [0.7 0.3 0.7] [0.3 0.7 0.3] [0.7 0.3 0.7]};
        
        % hgp
        %     figure()
        figure('Visible','Off');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 5])
        
        hold on
        
        ut_hand(1) = subplot(1,2,1);
        for iev = 3:4
            
            hold on
            
            mdat = squeeze(ft_dat_hgp(nchan_hgp(ic),:,ev_ind_hgp{iev}))';
            ev_mean = mean(mdat);
            
            plot(epoch_data_hgp_preproc.time{1},ev_mean,'Color',col{iev},'LineWidth',2);
            
%             % Deviation
            lowb  = ev_mean-(std(mdat)/sqrt(size(mdat,1)));
            highb = ev_mean+(std(mdat)/sqrt(size(mdat,1)));
%             
%             % Add measures of variance
%             patch([epoch_data_hgp_preproc.time{1} epoch_data_hgp_preproc.time{1}(end:-1:1)],[highb lowb(end:-1:1)],col_dev{iev})
%             chldrn = get(gca,'Children');
%             lnchldrn = findobj(gca,'Type','line');
%             [~,iH] = sort(ismember(chldrn,lnchldrn),1,'descend');
%             set(gca,'Children',chldrn(iH));
            
            % For plotting purposes to show the highest value
            maxy(1,iev) = max(highb);
            miny(1,iev) = min(lowb);
        end
        
        % Stats
        maxy_curr = max(maxy(:));
        miny_curr = min(miny(:));
        
        stats = stat_data_sa_semantic.stats(1).prob(nchan_hgp(ic),:) < 0.05;
        
        if stats(1) == 1; stats(1) = 0; end
        if stats(end) == 1; stats(end) = 0; end
        sig_periods = ts_crossing(stats(:),[],0.5); % - EJK change to 'crossing.m'
        sig_1st = sig_periods(1:2:end);
        sig_2nd = sig_periods(2:2:end);
        
        lnchldrn = findobj(gca,'Type','line');
        patchchldrn = findobj(gca,'Type','patch');
        combchldrn = [lnchldrn ; patchchldrn];
        
        for isig = 1:length(sig_1st)
            ftime = epoch_data_hgp_preproc.time{1}(sig_1st(isig):sig_2nd(isig)) ; btime = epoch_data_hgp_preproc.time{1}(sig_2nd(isig):-1:sig_1st(isig));
            patch([ftime btime],[repmat(maxy_curr,1,sig_2nd(isig)-sig_1st(isig)+1) repmat(miny_curr,1,sig_2nd(isig)-sig_1st(isig)+1)],[0.7 0.7 0.7])
            chldrn = get(gca,'Children');
            [~,iH] = sort(ismember(chldrn,combchldrn),1,'descend');
            set(gca,'Children',chldrn(iH));
        end
        
        % Remove Axis Information
        set(gca,'YTickLabel',[])
        set(gca,'XTickLabel',[])
        
        hold off
        
        ut_hand(2) = subplot(1,2,2);
        for iev = 1:2
            
            hold on
            
            mdat = squeeze(ft_dat_hgp(nchan_hgp(ic),:,ev_ind_hgp{iev}))';
            ev_mean = mean(mdat);
            
            plot(epoch_data_hgp_preproc.time{1},ev_mean,'Color',col{iev},'LineWidth',2);
%             
%             % Deviation
            lowb  = ev_mean-(std(mdat)/sqrt(size(mdat,1)));
            highb = ev_mean+(std(mdat)/sqrt(size(mdat,1)));
%             
%             % Add measures of variance
%             patch([epoch_data_hgp_preproc.time{1} epoch_data_hgp_preproc.time{1}(end:-1:1)],[highb lowb(end:-1:1)],col_dev{iev})
%             chldrn = get(gca,'Children');
%             lnchldrn = findobj(gca,'Type','line');
%             [~,iH] = sort(ismember(chldrn,lnchldrn),1,'descend');
%             set(gca,'Children',chldrn(iH));
            
            % For plotting purposes to show the highest value
            maxy(1,iev) = max(highb);
            miny(1,iev) = min(lowb);
        end
        
        % Stats
        maxy_curr = max(maxy(:));
        miny_curr = min(miny(:));
        
        stats = stat_data_sz_semantic.stats(1).prob(nchan_hgp(ic),:) < 0.05;
        
        if stats(1) == 1; stats(1) = 0; end
        if stats(end) == 1; stats(end) = 0; end
        sig_periods = ts_crossing(stats(:),[],0.5); % - EJK change to 'crossing.m'
        sig_1st = sig_periods(1:2:end);
        sig_2nd = sig_periods(2:2:end);
        
        lnchldrn = findobj(gca,'Type','line');
        patchchldrn = findobj(gca,'Type','patch');
        combchldrn = [lnchldrn ; patchchldrn];
        
        for isig = 1:length(sig_1st)
            ftime = epoch_data_hgp_preproc.time{1}(sig_1st(isig):sig_2nd(isig)) ; btime = epoch_data_hgp_preproc.time{1}(sig_2nd(isig):-1:sig_1st(isig));
            patch([ftime btime],[repmat(maxy_curr,1,sig_2nd(isig)-sig_1st(isig)+1) repmat(miny_curr,1,sig_2nd(isig)-sig_1st(isig)+1)],[0.7 0.7 0.7])
            chldrn = get(gca,'Children');
            [~,iH] = sort(ismember(chldrn,combchldrn),1,'descend');
            set(gca,'Children',chldrn(iH));
        end
        
        % Remove Axis Information
        set(gca,'YTickLabel',[])
        set(gca,'XTickLabel',[])
        
        linkaxes(ut_hand,'xy');
        xlim([epoch_data_hgp_preproc.time{1}(1) epoch_data_hgp_preproc.time{1}(end)])
        
        maxy = max(max(maxy));
        miny = min(min(miny));
        ylim([miny-miny*0.1 maxy+maxy*0.1])
        
        tick_y = (miny-miny*0.1):((abs(miny-miny*0.1)+abs(maxy+maxy*0.1))/3):(maxy+maxy*0.1);
        tick_x = epoch_data_hgp_preproc.time{1}(1):0.2:epoch_data_hgp_preproc.time{1}(end);
        
        set(ut_hand(:),'YTick',tick_y)
        set(ut_hand(:),'XTick',tick_x)
        
        linkaxes(ut_hand,'off');
        
        % Label first left axis
        set(ut_hand(1),'YTickLabel',tick_y,'FontSize',21)
        set(ut_hand(1),'XTickLabel',tick_x,'FontSize',21)
        
        % Add stimulus line
        for il = 1:length(ut_hand)
            subplot(1,2,il);
            line('XData',[0 0],'YData',[floor(miny-miny*0.1) ceil(maxy+maxy*0.1)],'LineWidth',2,'Color','k');
        end
        
         labe = strrep(epoch_data_lfp_preproc.label{nchan_lfp(ic)},'2>',''); 
        print(gcf,[indir '/' subj '/' 'vis' '/' subj '_' labe '_hgp_semantic'],'-dpng','-r500')
        close all
        
        cfg = [];
        cfg.return_events = 1;
        epoch_data_hgp_preproc = ft_redefine_events(cfg,epoch_data_hgp_preproc);
    
    end
    
end