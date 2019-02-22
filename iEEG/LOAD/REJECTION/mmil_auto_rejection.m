% cfg.chan     = 'all' (default) or index %% - EJK needs to be added
% cfg.thresh   = 'liberal' (0.99; default), 'conservative' (0.95), 'stringent' (0.90) (input single value between 0 and 1)
%     Set with number
% cfg.measures = 'all' (default), {'time' 'Hz2-12' 'Hz14-55' 'Hz70-170' 'variance'}
% cfg.dsamp    = apply downsampling to data
% cfg.epoch.[needed for epoching] = cfg needs to have all the necessary information to create epochs
% cfg.outdir   = directory to print figures to
% cfg.prefix   = name to prefix to the figures that are output
% cfg.pad      = numer of how much time on each side to trim (is s)
% cfg.plot     = plot the data? (1 is default)
% cfg.pct_chn  =  set number of channels which must share an artifact to remove from all channels (for use with meg settings)
% cfg.ext_chn  =  set extreme value to remove trial for all channels, regardless of how many channels have the artifact (for use with meg settings)

function ft_dat = auto_rej(fcfg,ft_dat)

if isfield(ft_dat,'trial');
    is_td = 1; is_fd = 0; is_pw = 0;
elseif isfield(ft_dat,'fourierspctrm');
    is_td = 0; is_fd = 1; is_pw = 0;
elseif isfield(ft_dat,'power');
    is_fd = 0; is_td = 0; is_pw = 1;
end

%% Set up Datatypes
if (~isfield(fcfg,'measures') || sum(strcmpi(fcfg.measures,'all'))) && is_td
    fcfg.measures = {'time' 'Hz2-12' 'Hz14-55' 'variance'};
elseif  (~isfield(fcfg,'measures') || sum(strcmpi(fcfg.measures,'all'))) && is_td
    fcfg.measures = {'Hz2-8' 'Hz9-12' 'Hz14-28' 'Hz30-55'};
end
if ~isfield(fcfg,'thresh');
    fcfg.thresh = repmat(0.99,1,numel(fcfg.measures));
elseif numel(fcfg.thresh) ~= numel(fcfg.measures);
    fcfg.thresh = repmat(fcfg.thresh,1,numel(fcfg.measures));
end
if isfield(fcfg,'pct_chn') && numel(fcfg.pct_chn) ~= numel(fcfg.measures);
    fcfg.pct_chn = repmat(fcfg.pct_chn,1,numel(fcfg.measures));
end
if isfield(fcfg,'ext_chn') && numel(fcfg.ext_chn) ~= numel(fcfg.measures);
    fcfg.ext_chn = repmat(fcfg.ext_chn,1,numel(fcfg.measures));
end

if ~isfield(fcfg,'chan') || strcmpi(fcfg.chan,'all'); fcfg.chan = 1:numel(ft_dat.label); end
if ~isfield(fcfg,'plot'); fcfg.plot = 1; end
if ~isdir(fcfg.outdir); mkdir(fcfg.outdir); end
if ~isfield(fcfg,'dsamp'); fcfg.dsamp = 0; end
if ~isfield(fcfg,'pad');   fcfg.pad   = 0; end

if is_td; smp_rte = ft_dat.fsample; elseif is_fd || is_pw; smp_rte = round(1 / (ft_dat.time(2) - ft_dat.time(1))); end

%% Load & Downsample
if fcfg.dsamp
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq   = 100;
    ft_dat      = ft_preprocessing(cfg,ft_dat);
    
    cfg = [];
    cfg.detrend = 'no';
    cfg.resamplefs = 250;
    ft_dat = ft_resampledata(cfg,ft_dat);
end

%% Create Waves
prc_dat = cell(1,numel(fcfg.measures));

if is_td
    for iMS = 1:numel(fcfg.measures)
        
        if strcmpi(fcfg.measures(iMS),'time') || strcmpi(fcfg.measures(iMS),'time-1')
            if isfield(fcfg,'epoch')
                cfg = [];
                cfg.lpfilter = 'yes';
                cfg.lpfreq   = [45];
                prc_dat{iMS} = ft_preprocessing(cfg,ft_dat);
            else
                prc_dat{iMS} = ft_dat;
            end
        elseif strcmpi(fcfg.measures(iMS),'variance')
            if isfield(fcfg,'epoch')
                cfg = [];
                cfg.lpfilter = 'yes';
                cfg.lpfreq   = [45];
                prc_dat{iMS} = ft_preprocessing(cfg,ft_dat);
            else
                prc_dat{iMS} = ft_dat;
            end
            
        else
            
            cfg = [];
            cfg.bpfilter = 'yes';
            
            if strcmpi(fcfg.measures(iMS),'Hz2-12')
                cfg.bpfreq   = [2 12];
            elseif strcmpi(fcfg.measures(iMS),'Hz14-55')
                cfg.bpfreq   = [14 55];
            end
            
            prc_dat{iMS} = ft_preprocessing(cfg,ft_dat);
            
            cfg = [];
            cfg.hilbert  = 'abs';
            prc_dat{iMS} = ft_preprocessing(cfg,prc_dat{iMS});
            
        end
        
    end
    
elseif is_fd
    
    prc_dat = cell(1,numel(fcfg.measures)*numel(fcfg.frq_rng));
    
    cnt = 1;
    
    for iFR = 1:numel(fcfg.frq_rng)
        for iMS = 1:numel(fcfg.measures)
            [~,beg_frq] = min(abs(ft_dat.freq-fcfg.frq_rng{iFR}(1))); [~,end_frq] = min(abs(ft_dat.freq-fcfg.frq_rng{iFR}(2)));
            
            prc_dat{cnt} = squeeze(mean(abs(ft_dat.fourierspctrm(:,:,beg_frq:end_frq,:)),3));
            
            if numel(size(prc_dat{cnt})) < 3; prc_dat{cnt} = reshape(prc_dat{cnt},[size(prc_dat{cnt},1) 1 size(prc_dat{cnt},2)]); end
            
            cnt = cnt + 1;
        end
    end
    
elseif is_pw
    
    prc_dat = cell(1,numel(fcfg.measures)*numel(fcfg.frq_rng));
    
    cnt = 1;
    
    for iFR = 1:numel(fcfg.frq_rng)
        for iMS = 1:numel(fcfg.measures)
            [~,beg_frq] = min(abs(ft_dat.freq-fcfg.frq_rng{iFR}(1))); [~,end_frq] = min(abs(ft_dat.freq-fcfg.frq_rng{iFR}(2)));
            
            prc_dat{cnt} = squeeze(mean(abs(ft_dat.power(:,:,beg_frq:end_frq,:)),3));
            
            cnt = cnt + 1;
        end
    end
    
    
end

%% Epoch Data
if isfield(fcfg,'epoch')
    
    cfg = [];
    epc_fld = fieldnames(fcfg.epoch);
    for iEP = 1:numel(epc_fld)
        cfg.(epc_fld{iEP}) = fcfg.epoch.(epc_fld{iEP});
    end
    cfg.fsample             = ft_dat.fsample;
    cfg                     = ft_definetrial(cfg);
    
    % segment data
    trl = cfg.trl;
    cfg = [];
    cfg.trl = trl;
    
    dat_hld = cell(1,numel(fcfg.measures));
    for iME = 1:numel(prc_dat)
        prc_dat{iME} = ft_redefinetrial(cfg,prc_dat{iME}); dat_hld{iME} =  permute(cat(3,prc_dat{iME}.trial{:}),[3 1 2]);
        poi_cut      = prc_dat{1}.fsample * fcfg.pad;
        dat_hld{iME}(:,:,poi_cut+1:end-poi_cut) = [];
    end
    clear prc_dat;
    
else
    
    dat_hld = cell(1,numel(prc_dat));
    for iME = 1:numel(prc_dat)
        if is_td; dat_hld{iME} = permute(cat(3,prc_dat{iME}.trial{:}),[3 1 2]); elseif is_fd; dat_hld{iME} = prc_dat{iME}; elseif is_pw; dat_hld{iME} = prc_dat{iME}; end
        poi_cut      = round(smp_rte * fcfg.pad);
        if poi_cut ~= 0; dat_hld{iME}(:,:,[1:poi_cut+1 end-poi_cut:end]) = []; end
    end
    
    if is_fd
        frq_rng_hld   = fcfg.frq_rng(ones(1,numel(fcfg.measures)),:);
        frq_rng_hld   = frq_rng_hld(:)';
        fcfg.measures = repmat(fcfg.measures,1,numel(fcfg.frq_rng));
        fcfg.thresh   = repmat(fcfg.thresh,1,numel(fcfg.frq_rng));
    end
    
end

%% Loop Over All Channels to find bad events on each channel
if is_td
    z_scr = zeros(numel(ft_dat.label),numel(fcfg.measures),numel(prc_dat{1}.trial));
    z_rmv = zeros(numel(ft_dat.label),numel(fcfg.measures),numel(prc_dat{1}.trial));
    z_cut = zeros(numel(ft_dat.label),numel(fcfg.measures));
    
    ft_dat.cfg.badtrl      = zeros(numel(fcfg.measures)+1,numel(ft_dat.label),numel(prc_dat{1}.trial));
    ft_dat.cfg.badtrl_type = [fcfg.measures 'All'];
elseif is_fd
    z_scr = zeros(numel(ft_dat.label),numel(fcfg.measures),size(prc_dat{1},1));
    z_rmv = zeros(numel(ft_dat.label),numel(fcfg.measures),size(prc_dat{1},1));
    z_cut = zeros(numel(ft_dat.label),numel(fcfg.measures));
    
    ft_dat.cfg.badtrl      = zeros(numel(fcfg.measures)+1,numel(ft_dat.label),size(prc_dat{1},1));
    ft_dat.cfg.badtrl_type = [fcfg.measures 'All'];
elseif is_pw
    z_scr = zeros(numel(ft_dat.label),numel(fcfg.measures),size(prc_dat{1},1));
    z_rmv = zeros(numel(ft_dat.label),numel(fcfg.measures),size(prc_dat{1},1));
    z_cut = zeros(numel(ft_dat.label),numel(fcfg.measures));
    
    ft_dat.cfg.badtrl      = zeros(numel(fcfg.measures)+1,numel(ft_dat.label),size(prc_dat{1},1));
    ft_dat.cfg.badtrl_type = [cellfun(@(x,y) strcat(num2str(x(1)),'-',num2str(x(2)),' ',y),frq_rng_hld,fcfg.measures,'uni',0) 'All'];
end

% Setup variables
z_scr = zeros(numel(fcfg.chan),numel(fcfg.measures),size(dat_hld{1},1));
z_rmv = zeros(numel(fcfg.chan),numel(fcfg.measures),size(dat_hld{1},1));

num_mea = numel(fcfg.measures);

% Run on every channel
for iCH = fcfg.chan % parfor iCH = fcfg.chan
    
    fprintf('Working on Channel %i\n',iCH)
    
    if numel(unique(dat_hld{1}(:,iCH,:))) ~= 1
    
    for iZS = 1:num_mea
        
        if strcmpi(fcfg.measures(iZS),'time')
            
            dat = squeeze(dat_hld{iZS}(:,iCH,:));
            
            z_men_dat    = nanmean(nanmean(dat));
            z_std_dat    = nanmean(nanstd(dat));
            z_scr(iCH,iZS,:) = max(((dat-z_men_dat)/z_std_dat)');
            
            [max_val max_val_ind] = max(squeeze(z_scr(iCH,iZS,:)));
            std_val               = nanstd(squeeze(z_scr(iCH,iZS,:)));
            std_med               = nanmedian(squeeze(z_scr(iCH,iZS,:)));
            
            while max_val > (std_med + std_val * 10)
                
                z_rmv_hld = z_rmv(iCH,iZS,:); z_rmv_hld(1,1,max_val_ind) = max_val;
                z_scr_hld = z_scr(iCH,iZS,:); z_scr_hld(1,1,max_val_ind) = nan;
                
                z_rmv(iCH,iZS,:) = z_rmv_hld;
                z_scr(iCH,iZS,:) = z_scr_hld;
                
                [max_val max_val_ind] = max(squeeze(z_scr(iCH,iZS,:)));
                std_val               = nanstd(squeeze(z_scr(iCH,iZS,:)));
                std_med               = nanmedian(squeeze(z_scr(iCH,iZS,:)));
            end
            
        elseif strcmpi(fcfg.measures(iZS),'time-1')
            
            dat = squeeze(dat_hld{iZS}(:,iCH,:)) * -1;
            
            z_men_dat    = nanmean(nanmean(dat));
            z_std_dat    = nanmean(nanstd(dat));
            z_scr(iCH,iZS,:) = max(((dat-z_men_dat)/z_std_dat)');
            
            [max_val,max_val_ind] = max(squeeze(z_scr(iCH,iZS,:)));
            std_val               = nanstd(squeeze(z_scr(iCH,iZS,:)));
            std_med               = nanmedian(squeeze(z_scr(iCH,iZS,:)));
            
            while max_val > (std_med + std_val * 10)
                
                z_rmv_hld = z_rmv(iCH,iZS,:); z_rmv_hld(1,1,max_val_ind) = max_val;
                z_scr_hld = z_scr(iCH,iZS,:); z_scr_hld(1,1,max_val_ind) = nan;
                
                z_rmv(iCH,iZS,:) = z_rmv_hld;
                z_scr(iCH,iZS,:) = z_scr_hld;
                
                [max_val,max_val_ind] = max(squeeze(z_scr(iCH,iZS,:)));
                std_val               = nanstd(squeeze(z_scr(iCH,iZS,:)));
                std_med               = nanmedian(squeeze(z_scr(iCH,iZS,:)));
            end
            
        elseif strcmpi(fcfg.measures(iZS),'variance')
            
            dat     = squeeze(dat_hld{iZS}(:,iCH,:));
            dat     = var(dat');
            
            z_scr(iCH,iZS,:) = dat;
            
            [max_val,max_val_ind] = max(squeeze(z_scr(iCH,iZS,:)));
            std_val               = nanstd(squeeze(z_scr(iCH,iZS,:)));
            std_med               = nanmedian(squeeze(z_scr(iCH,iZS,:)));
            
            while max_val > (std_med + std_val * 10)
                
                z_rmv_hld = z_rmv(iCH,iZS,:); z_rmv_hld(1,1,max_val_ind) = max_val;
                z_scr_hld = z_scr(iCH,iZS,:); z_scr_hld(1,1,max_val_ind) = nan;
                
                z_rmv(iCH,iZS,:) = z_rmv_hld;
                z_scr(iCH,iZS,:) = z_scr_hld;
                
                [max_val,max_val_ind] = max(squeeze(z_scr(iCH,iZS,:)));
                std_val               = nanstd(squeeze(z_scr(iCH,iZS,:)));
                std_med               = nanmedian(squeeze(z_scr(iCH,iZS,:)));
            end
            
        end
        
        if mean(squeeze(z_scr(iCH,iZS,:))) < 10^-4 && nanmean(squeeze(z_scr(iCH,iZS,:))) > -10^-4; z_scr(iCH,iZS,:) = z_scr(iCH,iZS,:) * 1/nanmean(squeeze(z_scr(iCH,iZS,:))) ; end
        
        % Calculate Cutoffs
        z_cut(iCH,iZS) = find_thresh(squeeze(z_scr(iCH,iZS,:)),fcfg.thresh(iZS));
        
        % Find Bad Indices
        rmv_nan     = z_scr(iCH,iZS,:);
        rmv_nan(find(isnan(rmv_nan))) = 100;
        z_scr(iCH,iZS,:) = rmv_nan;
        
        bad_ind_hld = squeeze(z_scr(iCH,iZS,:));
        bad_ind     = squeeze(z_scr(iCH,iZS,:)<z_cut(iCH,iZS));
        bad_ind_hld(bad_ind) = 0;
        
        z_rmv(iCH,iZS,:) = bad_ind_hld;
        
    end
    
    end
end

if isfield(fcfg,'pct_chn')
    for iMA = 1:numel(fcfg.measures)
        % Look for Channel Agreement
        [~,ic] = find(squeeze(z_rmv(:,iMA,:)));
        ic_pct = tabulate(ic); ic_pct = ic_pct(ic_pct(:,2)>size(ft_dat.cfg.badtrl,2)*fcfg.pct_chn(iMA),1);
        % Look for extreme Values
        ic_ext     = find(squeeze(z_scr(:,iMA,:))>fcfg.ext_chn(iMA));
        [~,ic_ext] = ind2sub(size(squeeze(z_scr(:,iMA,:))),ic_ext);
        ic_ext     = unique(ic_ext);
        
        ic_ovr{iMA} = unique([ic_pct' ic_ext']);
        % Apply to Data
        if ~isempty(ic_ovr{iMA});
            rmv_ind = sub2ind(size(ft_dat.cfg.badtrl),ones(1,numel(ft_dat.label)*numel(ic_ovr{iMA}))*iMA,repmat(1:numel(ft_dat.label),1,numel(ic_ovr{iMA})),cell2mat(arrayfun(@(x, y) repmat(x, [1 y]),ic_ovr{iMA},repmat(numel(ft_dat.label),1,numel(ic_ovr{iMA})),'UniformOutput',0)));
            ft_dat.cfg.badtrl(rmv_ind) = 1;
        end
        
    end
    ic_ovr = unique([ic_ovr{:}]);
    rmv_ind = sub2ind(size(ft_dat.cfg.badtrl),ones(1,numel(ft_dat.label)*numel(ic_ovr))*iMA+1,repmat(1:numel(ft_dat.label),1,numel(ic_ovr)),cell2mat(arrayfun(@(x, y) repmat(x, [1 y]),ic_ovr,repmat(numel(ft_dat.label),1,numel(ic_ovr)),'UniformOutput',0)));
    ft_dat.cfg.badtrl(rmv_ind) = 1;
elseif ~isfield(fcfg,'pct_chn')
    % Save Removed Trials
    for iCH = fcfg.chan
    ft_dat.cfg.badtrl(1:numel(fcfg.measures),iCH,:) = squeeze(z_rmv(iCH,:,:));
    ft_dat.cfg.badtrl(numel(fcfg.measures)+1,iCH,:) = [squeeze(sum(z_rmv(iCH,:,:)))>0]';
    end
end

%% Make Plots of Findings
if fcfg.plot
    % Individual Channel Plots
    for iCH = fcfg.chan % parfor iCH = fcfg.chan
        
        fprintf('Printing Channel %i\n',iCH)
        
        hFig = figure('Visible','off');
        set(hFig,'Units','Normalized','Position', [0 0 1 1])
        for iPL = 1:numel(fcfg.measures)
            subplot(ceil((numel(fcfg.measures)+1)/3),3,iPL); hold on;
            
            if ~strcmpi(fcfg.measures{iPL},'variance')
                if is_td || is_pw
                    acc_trl = find(~z_rmv(iCH,iPL,:)); rej_trl = find(z_rmv(iCH,iPL,:));
                    
                    for iLN = 1:numel(acc_trl); plot(squeeze(dat_hld{iPL}(acc_trl(iLN),iCH,:))); end;
                    for iLN = 1:numel(rej_trl); plot(squeeze(dat_hld{iPL}(rej_trl(iLN),iCH,:)),'r'); end;
                    for iLN = 1:numel(rej_trl); plot(squeeze(dat_hld{iPL}(rej_trl(iLN),iCH,:)),'r'); end;
                    
                    if ~is_pw
                        title([fcfg.measures{iPL} ' ' 'rej: ' num2str(numel(rej_trl))])
                    elseif is_pw
                        title([ft_dat.cfg.badtrl_type{iPL} ' ' 'rej: ' num2str(numel(rej_trl))])
                    end
                    
                    
                elseif is_fd
                    acc_trl     = find(~z_rmv(iCH,iPL,:)); rej_trl_kpt = find(z_rmv(iCH,iPL,:)); rej_trl = find(ft_dat.cfg.badtrl(end,iCH,:));
                    rej_trl_kpt = setdiff(rej_trl_kpt,rej_trl); rej_trl =  setdiff(rej_trl,acc_trl);
                    ovr_rej     = find(ft_dat.cfg.badtrl(end,iCH,:)); ovr_rej = setdiff(ovr_rej,rej_trl);
                    
                    for iLN = 1:numel(acc_trl); plot(squeeze(dat_hld{iPL}(acc_trl(iLN),iCH,:))); end;
                    for iLN = 1:numel(rej_trl_kpt); plot(squeeze(dat_hld{iPL}(rej_trl_kpt(iLN),iCH,:)),'c'); end;
                    for iLN = 1:numel(rej_trl); plot(squeeze(dat_hld{iPL}(rej_trl(iLN),iCH,:)),'r'); end;
                    for iLN = 1:numel(ovr_rej); plot(squeeze(dat_hld{iPL}(ovr_rej(iLN),iCH,:)),'m'); end;
                    
                    title([fcfg.measures{iPL} ' ' 'Rej: ' num2str(numel(rej_trl)) ' ' 'Kpt: ' num2str(numel(rej_trl_kpt)) ' ' 'OvrRej: ' num2str(numel(ovr_rej)) ' ' 'Acc: ' num2str(numel(acc_trl))]);
                    
                end
            else
                
                plt_var = squeeze(z_scr(iCH,iPL,:));
                
                try
                    [y,x] = hist(plt_var,300);
                    x1 = x(x<z_cut(iCH,iPL)); y1 = y(x<z_cut(iCH,iPL));
                    x2 = x(x>z_cut(iCH,iPL)); y2 = y(x>z_cut(iCH,iPL));
                    bar(x1,y1,'b'); hold on; bar(x2,y2,'r');
                    xlimit = get(gca,'xlim'); set(gca,'xlim',[0 xlimit(2)]);
                    vline(z_cut(iCH,iPL),'r');
                    
                catch err
                end

                if ~is_pw
                    title([fcfg.measures{iPL} ' ' 'rej: ' num2str(sum(y2))])
                elseif is_pw
                    title([ft_dat.cfg.badtrl_type{iPL} ' ' 'rej: ' num2str(sum(y2))])
                end
            end
            
        end
        
        % Overal Rejection Plot
        subplot(ceil((numel(fcfg.measures)+1)/3),3,iPL+1); hold on;
        
        acc_trl = find(~ft_dat.cfg.badtrl(end,iCH,:)); rej_trl_kpt = unique(find(z_rmv(iCH,:,:))); [~,rej_trl_kpt] = ind2sub(size(squeeze(z_rmv(iCH,:,:))),rej_trl_kpt); rej_trl = find(ft_dat.cfg.badtrl(end,iCH,:));
        rej_trl_kpt = setdiff(rej_trl_kpt,rej_trl); acc_trl =  setdiff(acc_trl,rej_trl_kpt);
        
        if is_td || is_pw
            
            for iLN = 1:numel(acc_trl); plot(squeeze(dat_hld{iPL}(acc_trl(iLN),iCH,:))); end;
            for iLN = 1:numel(rej_trl_kpt); plot(squeeze(dat_hld{iPL}(rej_trl_kpt(iLN),iCH,:)),'c'); end;
            for iLN = 1:numel(rej_trl); plot(squeeze(dat_hld{iPL}(rej_trl(iLN),iCH,:)),'r'); end;
            
        end
        
        title(['Rej: ' num2str(numel(rej_trl)) ' ' 'Kpt: ' num2str(numel(rej_trl_kpt)) ' ' 'Acc: ' num2str(numel(acc_trl))]);
        
        print(gcf,[fcfg.outdir '/' fcfg.prefix '_' ft_dat.label{iCH} '.png'],'-dpng')
        close all
        
    end
    
    
end

end

%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
