function [stt_dat] = mmil_ieeg_statistics5(cfg,dat)
%% Notes
% WILL ADD MONTE CARLO
% WILL ADD TimeFrequency Stats
% CAN ADD 2x2 ANOVAS
% CAN ADD Time Window to perform stats over cfg.stt_tme
% CAN ADD fdr correction for timewindow statistics

%% BEGIN
if ~isfield(cfg,'eve');     error('Must specify eve'); end

if ~isfield(cfg,'alt_eve'); error('Must specify alt_eve'); end
if ~isfield(cfg,'cor_mth'); cfg.cor_mth = ''; end
if ~isfield(cfg,'alp'); cfg.alp = .05; end
if ~isfield(cfg,'trl'); cfg.trl = 1:size(cfg.dta,3); end
if ~isfield(cfg,'chn'); cfg.chn = 1:numel(cfg.lbl); end
if ~isfield(cfg,'cfg.num_rsm'); cfg.num_rsm = 1000; end

if ~isfield(cfg,'fdr_typ'); cfg.fdr_typ = 'tme'; end

%% One-Way ANOVA (Unbalanced Design) and FUTURE Two-Way ANOVA
% Trials of interest
if strcmpi(cfg.typ,'trial')
    
    if ~isfield(cfg,'tme_win') && (strcmpi(cfg.cor_mth,'fdr') || strcmpi(cfg.cor_mth,'mnt_car') || strcmpi(cfg.cor_mth,''))  % TIMEPOINT-by-TIMEPOINT
        pvl_mtx = zeros([numel(cfg.chn) numel(cfg.tme{1})]);  % Creating an empty matrix to store p-values %%%%%%%%%%%%%%%%%%%
        fst_mtx = zeros([numel(cfg.chn) numel(cfg.tme{1})]);
    end
    
    if isfield(cfg,'tme_win')
        pvl_mtx = zeros([numel(cfg.chn) 1]);
        fst_mtx = zeros([numel(cfg.chn) 1]);
    end
    
    if strcmpi(cfg.cor_mth,'mnt_car') || numel(cfg.eve)==1
        pvl_rsm_mtx = zeros([numel(cfg.chn) numel(cfg.tme{1}) cfg.num_rsm]);
        fst_rsm_mtx = zeros([numel(cfg.chn) numel(cfg.tme{1}) cfg.num_rsm]);
    end
    
    for iC = 1:numel(cfg.chn)
        
        % remove NAN's from cfg.trl
        cfg.trl = intersect(cfg.trl,find(~isnan(squeeze(cfg.dta(iC,1,:)))));
        
        % events
        eve_cat = [];
        eve_ind = [];
        for iE = 1:numel(cfg.eve)
            try eve_ind = [eve_ind intersect(find(cfg.eve_hld==cfg.eve(iE)),cfg.trl)']; catch; eve_ind = [eve_ind ; intersect(find(cfg.eve_hld==cfg.eve(iE)),cfg.trl)']; end
            eve_cat = [eve_cat ones(1,numel(intersect(find(cfg.eve_hld==cfg.eve(iE)),cfg.trl)))*iE];
        end
        dat_hld = cfg.dta(iC,:,eve_ind);
        
        if isfield(cfg,'tme_win')
            if numel(cfg.tme_win) == 2
                cfg.tme_win{1} = dsearchn(cfg.tme{1}',cfg.tme_win{1}')'; cfg.tme_win{2} = dsearchn(cfg.tme{2}',cfg.tme_win{2}')';
                dat_hld = [squeeze(nanmean(dat_hld(:,cfg.tme_win{1}(1):cfg.tme_win{1}(2),:),2)) ; squeeze(nanmean(dat_hld(:,cfg.tme_win{2}(1):cfg.tme_win{2}(2),:),2))]';
                eve_cat = [eve_cat eve_cat*2];
            else
                cfg.tme_win{1} = dsearchn(cfg.tme{1}',cfg.tme_win{1}')';
                dat_hld = squeeze(nanmean(dat_hld(:,cfg.tme_win{1}(1):cfg.tme_win{1}(2),:),2))';
            end
        end
        
        % Run Stats
        if ~isfield(cfg,'tme_win') && (strcmpi(cfg.cor_mth,'fdr') || strcmpi(cfg.cor_mth,'mnt_car') || strcmpi(cfg.cor_mth,''))  % TIMEPOINT-by-TIMEPOINT
            
            if strcmpi(cfg.cor_mth,'mnt_car') || numel(cfg.eve)==1 % Set up shuffle if monte carlo %%%%%%%%%%%%%%%%%%%
                scfg = [];
                scfg.resampling       = 'permutation';
                scfg.numrandomization = cfg.num_rsm;
                dat_eve_rsm_hld       = resampledesign(scfg,eve_cat);
                for iS = 1:cfg.num_rsm; dat_eve_rsm_hld(iS,:) = eve_cat(dat_eve_rsm_hld(iS,:)); end
            end
            
            if numel(cfg.eve) == 1  % TTEST OVERALL ACTIVATON %%%%%%%%%%%%%%%%%%%
                
                statfun = @statfun_indepsamplesT; % OBSERVATION STATS %%%%%%%%%%%%%%%%%%%
                scfg = [];
                scfg.ivar = 1;
                scfg.tail = 0;
                scfg.computeprob = 'yes';
                obs  = statfun(scfg,squeeze(dat_hld(iC,:,:)),eve_cat);
                
                pvl_mtx(iC,:) = obs.prob;
                fst_mtx(iC,:) = obs.stat;
                
                if strcmpi(cfg.cor_mth,'mnt_car') || numel(cfg.eve)==1 % SHUFFLE STATS %%%%%%%%%%%%%%%%%%%
                    
                    cls_crt_val = icdf('T',[cfg.alp/2 (1-cfg.alp/2)],size(dat_hld,3)-1);
                    
                    for iS = 1:cfg.num_rsm;
                        statfun = @statfun_indepsamplesT;
                        scfg = [];
                        scfg.ivar = 1;
                        scfg.tail = 0;
                        scfg.computeprob = 'yes';
                        obs  = statfun(scfg,squeeze(dat_hld(iC,:,:)),dat_eve_rsm_hld(iS,:));
                        
                        pvl_rsm_mtx(iC,:,iS) = obs.prob;
                        fst_rsm_mtx(iC,:,iS) = obs.stat;
                    end
                end
                
            elseif numel(cfg.eve) > 1 && ~isfield(cfg,'tme_win') % ANOVA DIFFERENCES IN ACTIVATON %%%%%%%%%%%%%%%%%%%
                
                if numel(unique(eve_cat)) == numel(cfg.eve)
                    statfun = @statfun_indepsamplesF;
                    scfg = [];
                    scfg.ivar = 1;
                    scfg.computeprob = 'yes';
                    obs  = statfun(scfg,squeeze(dat_hld(iC,:,:)),eve_cat);
                    
                    pvl_mtx(iC,:) = obs.prob;
                    fst_mtx(iC,:) = obs.stat;
                    
                    if strcmpi(cfg.cor_mth,'mnt_car') % SHUFFLE STATS %%%%%%%%%%%%%%%%%%%
                        
                        cls_crt_val = icdf('F',1-cfg.alp,numel(cfg.eve)-1,size(dat_hld,3)-1);
                        tail        = 1;
                        
                        for iS = 1:cfg.num_rsm;
                            statfun = @statfun_indepsamplesF;
                            scfg = [];
                            scfg.ivar = 1;
                            scfg.computeprob = 'yes';
                            obs  = statfun(scfg,squeeze(dat_hld(iC,:,:)),dat_eve_rsm_hld(iS,:));
                            
                            pvl_rsm_mtx(iC,:,iS) = obs.prob;
                            fst_rsm_mtx(iC,:,iS) = obs.stat;
                        end
                    end
                else
                    
                    pvl_mtx(iC,:) = ones(1,size(dat_hld,2));
                    fst_mtx(iC,:) = zeros(1,size(dat_hld,2));
                    
                    pvl_rsm_mtx(iC,:,:) = ones(1,size(dat_hld,2),cfg.num_rsm);
                    fst_rsm_mtx(iC,:,:) = zeros(1,size(dat_hld,2),cfg.num_rsm);
                    
                end
                
            elseif isfield(cfg,'tme_win') % TIMEWINDOW-by-TIMEWINDOW  %%%%%%%%%%%%%%%%%%%
                
                if numel(cfg.tme_win)==2 % OVERALL ACTIVATON %%%%%%%%%%%%%%%%%%%
                    
                    for iC = 1:numel(cfg.chn)
                        statfun = @statfun_indepsamplesT; % OBSERVATION STATS %%%%%%%%%%%%%%%%%%%
                        scfg = [];
                        scfg.ivar = 1;
                        scfg.tail = 0;
                        scfg.computeprob = 'yes';
                        obs  = statfun(scfg,squeeze(dat_hld(iC,:)),eve_cat);
                        
                        pvl_mtx(iC,:) = obs.prob;
                        fst_mtx(iC,:) = obs.stat;
                        
                        cfg.tme_win{1} = cfg.tme_win{2};
                    end
                    
                elseif numel(cfg.tme_win)==1 % DIFFERENCES IN ACTIVATON %%%%%%%%%%%%%%%%%%%
                    
                    statfun = @statfun_indepsamplesT; % OBSERVATION STATS %%%%%%%%%%%%%%%%%%%
                    scfg = [];
                    scfg.ivar = 1;
                    scfg.tail = 0;
                    scfg.computeprob = 'yes';
                    obs  = statfun(scfg,squeeze(dat_hld(iC,:)),eve_cat);
                    
                    pvl_mtx(iC,:) = obs.prob;
                    fst_mtx(iC,:) = obs.stat;
                    
                end
                
            end
            
        end
        
    end
    
elseif strcmpi(cfg.typ,'fourierspctrm')
    
    % Data Holder
    dat_hld = dat.fourierspctrm;
    
end

fprintf('Finished running Stats\n\n')

%% Correct P-Values (currently: fdr, shuffle, none)
if strcmpi(cfg.cor_mth,'fdr') || strcmpi(cfg.cor_mth,'mnt_car')
    
    if ~isfield(cfg,'tme_win') % TIMEPOINT-by-TIMEPOINT %%%%%%%%%%%%%%%%%%%
        
        if (strcmpi(cfg.cor_mth,'mnt_car') || strcmpi(cfg.cor_mth,'fdr')) && strcmpi(cfg.typ,'trial') % FDR TIMEDOMAIN %%%%%%%%%%%%%%%%%%%
            
            if strcmpi(cfg.cor_mth,'fdr')
                switch cfg.fdr_typ
                    
                    case 'chn_tme'
                        [~,fdr_msk] = fdr_cor([pvl_mtx(:)]',cfg.alp);
                        fdr_msk = vec2mat(fdr_msk,size(pvl_mtx,1))';
                        
                        if strcmpi(cfg.cor_mth,'mnt_car') && strcmpi(cfg.typ,'trial')
                            error('Not implemented')
                        end
                        
                    case 'chn'
                        fdr_msk=zeros(size(pvl_mtx))';
                        for iT = 1:size(pvl_mtx,2)
                            [tmp, fdr_msk(iT,:)]=fdr_cor([pvl_mtx(:,iT)]',cfg.alp);
                        end
                        fdr_msk = fdr_msk';
                        
                        if strcmpi(cfg.cor_mth,'mnt_car') && strcmpi(cfg.typ,'trial')
                            error('Not implemented')
                        end
                        
                    case 'tme'
                        fdr_msk=zeros(size(pvl_mtx));
                        for iC = 1:size(pvl_mtx,1)
                            [tmp, fdr_msk(iC,:)]=fdr_cor(pvl_mtx(iC,:),cfg.alp);
                        end
                        
                        if (strcmpi(cfg.cor_mth,'mnt_car') || numel(cfg.eve)==1) && strcmpi(cfg.typ,'trial')
                            fdr_rsm_msk=zeros(size(pvl_rsm_mtx));
                            for iC = 1:size(pvl_rsm_mtx,1)
                                for iS = 1:cfg.num_rsm;
                                    [tmp, fdr_rsm_msk(iC,:,iS)]=fdr_cor(squeeze(pvl_rsm_mtx(iC,:,iS)),cfg.alp);
                                end
                            end
                            pvl_rsm_msk(~fdr_rsm_msk) = 1;
                            fst_rsm_mtx(~fdr_rsm_msk) = 0;
                        end
                end
                
                pvl_mtx(~fdr_msk) = 1;
                fst_mtx(~fdr_msk) = 0;
            end
            
            if strcmpi(cfg.cor_mth,'mnt_car') || numel(cfg.eve)==1 % SHUFFLE TIMEDOMAIN %%%%%%%%%%%%%%%%%%%
                
                issource = 0;
                
                cls_crt_val = icdf('T',[cfg.alp/2 (1-cfg.alp/2)],size(dat_hld,3)-1);
                tail        = 0;
                
                for iC = 1:size(fst_mtx,1)
                    stt_rnd  = squeeze(fst_rsm_mtx(iC,:,:));
                    stt_obs  = squeeze(fst_mtx(iC,:))';
                    
                    scfg = [];
                    scfg.neighbors = [];
                    scfg.numrandomization = cfg.num_rsm;
                    scfg.channel          = cfg.lbl{iC};
                    scfg.clusterstatistic = 'maxsum';
                    scfg.dim              = [1 1 numel(cfg.tme{1})];
                    scfg.clusterthreshold = 'parametric';
                    scfg.clusteralpha     = 0.05;
                    scfg.clustercritval   = cls_crt_val;
                    scfg.tail             = tail;
                    scfg.clustertail      = tail;
                    scfg.feedback         = 'no';
                    [stat, scfg] = clusterstats(scfg, stt_rnd, stt_obs,'issource',issource);
                    
                    pvl_mtx(iC,:) = stat.prob; pvl_mtx(iC,pvl_mtx(iC,:)>cfg.alp) = 1;
                    
                end
                
            end
            
        elseif strcmpi(cfg.cor_mth,'fdr') && strcmpi(cfg.typ,'fourierspctrm') % FDR FREQUENCY %%%%%%%%%%%%%%%%%%%
            
            fprintf('Beginning fdr correction\n\n')
            
            switch cfg.fdr_typ
                
                case 'chn_tme'
                    [tmp,fdr_msk] = fdr_cor([pvl_mtx(:,iC,:)]',cfg.alp);
                    fdr_msk = vec2mat(fdr_msk,size(pvl_mtx,1))';
                case 'chn'
                    fdr_msk=zeros(size(pvl_mtx))';
                    for iT = 1:size(pvl_mtx,3)
                        [tmp, fdr_msk(:,:,iT)]=fdr_cor([pvl_mtx(:,:,iT)]',cfg.alp);
                    end
                    fdr_msk = fdr_msk';
                case 'tme'
                    fdr_msk=zeros(size(pvl_mtx));
                    for iC = 1:size(pvl_mtx,2)
                        [tmp, fdr_msk(:,iC,:)]=fdr_cor(pvl_mtx(:,iC,:),cfg.alp);
                    end
            end
            
            pvl_mtx(~fdr_msk) = 1;
            
        end
        
    elseif isfield(cfg,'tme_win') % TIMEWINDOW %%%%%%%%%%%%%%%%%%%
        
        pvl_mtx(~fdr_msk) = 1;
        
        pvl_mtx_hld = ones([numel(cfg.chn) numel(cfg.tme{1})]);
        pvl_mtx_hld(:,cfg.tme_win{1}(1):cfg.tme_win{1}(2)) = repmat(pvl_mtx,1,numel(cfg.tme_win{1}(1):cfg.tme_win{1}(2)));
        pvl_mtx = pvl_mtx_hld;
        
    end
    
elseif strcmpi(cfg.cor_mth,'')
    fprintf('No Corrections Specified\n\n')
    
    if ~isfield(cfg,'tme_win') % TIMEPOINT-by-TIMEPOINT
        
        pvl_mtx(pvl_mtx>cfg.alp) = 1;
        
    elseif isfield(cfg,'tme_win') % TIMEWINDOW
        
        pvl_mtx(pvl_mtx>cfg.alp) = 1;
        
        pvl_mtx_hld = ones([numel(cfg.chn) numel(cfg.tme{1})]);
        pvl_mtx_hld(:,cfg.tme_win{1}(1):cfg.tme_win{1}(2)) = repmat(pvl_mtx,1,numel(cfg.tme_win{1}(1):cfg.tme_win{1}(2)));
        pvl_mtx = pvl_mtx_hld;
        
    end
    
end

stt_dat = pvl_mtx;

fprintf('Finished corrections\n\n')

end
%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDR function code: created by Brian Quinn
function [fdrthresh,fdrmask] = fdr_cor(pval_vec,alpha)

pval_sorted = sort(pval_vec);
fdrpvec=((1:length(pval_vec))./length(pval_vec))*alpha;
fdrtest = pval_sorted-fdrpvec;
fdrp = find(fdrtest < 0,1,'last');
fdrmask = zeros(1,length(pval_vec));
if ~isempty(fdrp),
    fdrthresh = pval_sorted(fdrp);
    fdrmask(pval_vec < fdrthresh)=1;
else
    fdrthresh=0;
end

end
