function dat = ft_remove_common_noise2(cfg,dat)

if strcmpi(cfg.chn_nse_grp,'split');
    
    cfg = rmfield(cfg,'chn_nse_grp');
    
    cfg.chn_nse_grp_nme{1} = 'ecog';
    cfg.chn_nse_grp_nme{2} = 'noisy_ecog';
    cfg.chn_nse_grp_nme{3} = 'depth';
    cfg.chn_nse_grp_nme{4} = 'noisy_depth';
    
    ele_idn = mmil_readtext([cfg.clr_fld '/' 'ele_idn' '/' cfg.sbj_nme '/' cfg.dat_nme '_tmp.csv']);
    nse_chn =  mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'nse_chn'); nse_chn=nse_chn{1};
    
    cfg.chn_nse_grp{1} = find(cell2mat(ele_idn(:,3))==1);
    [cfg.chn_nse_grp{2},rmv_ind] = intersect(cfg.chn_nse_grp{1},nse_chn);
    cfg.chn_nse_grp{1}(rmv_ind) = [];
    
    cfg.chn_nse_grp{3} = find(cell2mat(ele_idn(:,3))==2);
    [cfg.chn_nse_grp{4},rmv_ind] = intersect(cfg.chn_nse_grp{3},nse_chn);
    cfg.chn_nse_grp{3}(rmv_ind) = [];
    
    cfg.spl_chn = cfg.chn_nse_grp;
    cfg.spl_nme = cfg.chn_nse_grp_nme;
    
elseif strcmpi(cfg.chn_nse_grp,'pedot');
    
    cfg = rmfield(cfg,'chn_nse_grp');
    
    ele_idn = mmil_readtext([cfg.clr_fld '/' 'ele_idn' '/' cfg.sbj_nme '/' cfg.dat_nme '.csv']);
    ele_gdd = mmil_readtext([cfg.clr_fld '/' 'ele_idn' '/' cfg.sbj_nme '/' cfg.dat_nme '_gdd_chn.csv']);
    
    cfg.chn_nse_grp{1} = find(cell2mat(ele_idn(:,3))==1 & cell2mat(ele_gdd(:,3))==1);
    cfg.chn_nse_grp{2} = find(cell2mat(ele_idn(:,3))==1 & cell2mat(ele_gdd(:,3))==0);
    
    cfg.spl_chn = cfg.chn_nse_grp;
    cfg.spl_nme = cfg.pre_fix;
    
    cfg.chn_nse_grp_nme = cfg.pre_fix;
    
else
    
    cfg.spl_chn         = cfg.chn_nse_grp;
    cfg.chn_nse_grp_nme = mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'cmn_nme');
    
end

dat.cfg.rmv_nse.spl_chn = cfg.spl_chn;
dat.cfg.rmv_nse.spl_nme = cfg.chn_nse_grp_nme;
dat.cfg.rmv_nse.pre_fix = cfg.pre_fix;

if ~isfield(cfg,'rnd_tme')
    
    trl_qtl     = [1:round(quantile(1:numel(dat.trial),0.1)):numel(dat.trial) numel(dat.trial)];
    if numel(trl_qtl)==10; trl_qtl     = [1:round(quantile(1:numel(dat.trial),0.09)):numel(dat.trial) numel(dat.trial)]; end
    for iQ = 1:numel(trl_qtl)-1
        rnd_trl{iQ,1} = randsample(trl_qtl(iQ):trl_qtl(iQ+1),10);
    end
    
    tme_ind = 1:numel(dat.time{1})-round(dat.fsample)-1;
    if isempty(tme_ind); 
        cfg.rnd_tme=ones(10,1); 
    else
        cfg.rnd_tme = sort(randsample(dat.time{1}(tme_ind),10));
    cfg.rnd_tme = dsearchn(dat.time{1}',cfg.rnd_tme');
    end
    
end

iGR = 1;
while iGR <= numel(cfg.chn_nse_grp)
    
    if ~isempty(cfg.chn_nse_grp{iGR})
        
        rmv_hld_ind = numel(cfg.chn_nse_grp) + 1;
        rmv_hld_nme = strsplit(cfg.pre_fix{iGR},cfg.dat_nme); rmv_hld_nme = [rmv_hld_nme{1} '_2'];
        
        for iQ = 1:numel(rnd_trl)
            for iT = 1:numel(rnd_trl{iQ})
                std_hld(iT) = nanstd(dat.trial{rnd_trl{iQ}(iT)}(cfg.chn_nse_grp{iGR},cfg.rnd_tme(iT)));
            end
            zzz = find(max(std_hld)==std_hld); if numel(zzz)>1; zzz = zzz(1); end
            cfg.rnd_trl(iQ) = rnd_trl{iQ}(zzz);
        end
        
        if ~isdir([cfg.out_dir '/' cfg.sbj_nme '/' cfg.pre_fix{iGR}]); mkdir([cfg.out_dir '/' cfg.sbj_nme '/' cfg.pre_fix{iGR}]); end
        
        ttt_dat = cat(3,dat.trial{:});
        ttt_men = squeeze(mean(ttt_dat(cfg.chn_nse_grp{iGR},:,:),1))';
        
        % Initial Plot
        col_plt = distinguishable_colors(numel(cfg.chn_nse_grp{iGR}));
        for iRN = 1:numel(cfg.rnd_tme)
            h(iRN) = figure('Visible','off'); subplot(3,1,1); hold on;
            for iP = 1:numel(dat.label(cfg.chn_nse_grp{iGR})); plot(dat.time{1},dat.trial{cfg.rnd_trl(iRN)}(cfg.chn_nse_grp{iGR}(iP),:),'Color',col_plt(iP,:)); end
            plot(dat.time{1},ttt_men(cfg.rnd_trl(iRN),:),'k','LineWidth',1.25);
            xlim([dat.time{1}(cfg.rnd_tme(iRN)) dat.time{1}(cfg.rnd_tme(iRN)+round(dat.fsample)-1)]);
        end
        for iRN = 1:numel(cfg.rnd_tme)
            set(0,'currentfigure',h(iRN)); subplot(3,1,2); hold on;
            plot(dat.time{1},ttt_men(cfg.rnd_trl(iRN),:),'k');
            xlim([dat.time{1}(cfg.rnd_tme(iRN)) dat.time{1}(cfg.rnd_tme(iRN)+round(dat.fsample)-1)]);
        end
        
        % If desired, remove rogue channels from plot
        if isfield(cfg,'rmv_chn') && cfg.rmv_chn == 1
            
            sat = 0;
            rmv_chn{iGR} = {};
            
            if exist([cfg.out_dir '/' cfg.sbj_nme '/' cfg.pre_fix{iGR} '/'  'rmv_num.mat'])
                load([cfg.out_dir '/' cfg.sbj_nme '/' cfg.pre_fix{iGR} '/'  'rmv_num.mat'])
                rmv_chn{iGR} = sve_num;
                [~,ind,~] = intersect(dat.label(cfg.chn_nse_grp{iGR}),sve_num);
                cfg.chn_nse_grp{iGR}(ind) = [];
                col_plt(ind,:) = [];
                ttt_men = squeeze(mean(ttt_dat(cfg.chn_nse_grp{iGR},:,:),1))';
            end
            
            if cfg.rmv_chn
                while sat == 0
                    
                    ddd = figure();
                    col_plt = distinguishable_colors(numel(cfg.chn_nse_grp{iGR}));
                    for iRN = 1:numel(cfg.rnd_tme)-1
                        subplot(5,2,iRN); hold off;
                        for iP = 1:numel(dat.label(cfg.chn_nse_grp{iGR})); plot(dat.time{1}+10e6*iRN,dat.trial{cfg.rnd_trl(iRN)}(cfg.chn_nse_grp{iGR}(iP),:),'Color',col_plt(iP,:));  hold on; end;
                        plot(dat.time{1}+10e6*iRN,ttt_men(cfg.rnd_trl(iRN),:),'k','LineWidth',1.25);
                        xlim([dat.time{1}(cfg.rnd_tme(iRN))+10e6*iRN dat.time{1}(cfg.rnd_tme(iRN)+round(dat.fsample)-2/2)+10e6*iRN]);
                        set(gca,'xticklabel',[])
                    end
                    subplot(5,2,10); plot((1:10)+10e6*(iRN+1),repmat(5,1,10)); hold on; line([5+10e6*(iRN+1) 5+10e6*(iRN+1)],[0 5]); text(4.5+10e6*(iRN+1),7,'Save','FontSize',20); text(2.5+10e6*(iRN+1),3,'Trash','FontSize',20); text(6.5+10e6*(iRN+1),3,'Trash All','FontSize',20); ylim([0 10]); set(gca,'xticklabel',[])
                    set(gcf,'Position',get(0,'ScreenSize'));
                    tightfig();
                    title([mmil_spec_char(cfg.sbj_nme,{'_'}) ' - ' mmil_spec_char(cfg.pre_fix{iGR},{'_'}) '  ' 'REMOVED: ' num2str(numel(rmv_chn{iGR}))],'FontSize',20)
                    
                    [plt_rmv,thr_rmv] = ginput(1); tme_rmv = plt_rmv-10e6*round(plt_rmv/10e6); plt_rmv = round(plt_rmv/10e6); if plt_rmv ~= 10; tme_rmv = dsearchn(dat.time{1}',tme_rmv); end;
                    %                     figure(); hist(dat.trial{cfg.rnd_trl(plt_rmv)}(cfg.chn_nse_grp{iGR},tme_rmv),100); vline(thr_rmv,'r')
                    %                     dat.time{1}(tme_rmv)
                    
                    if plt_rmv ~= 10;
                        
                        if thr_rmv < mean(dat.trial{cfg.rnd_trl(plt_rmv)}(cfg.chn_nse_grp{iGR},tme_rmv))
                            [row,col] = find(dat.trial{cfg.rnd_trl(plt_rmv)}(cfg.chn_nse_grp{iGR},tme_rmv) < thr_rmv);
                        else
                            [row,col] = find(dat.trial{cfg.rnd_trl(plt_rmv)}(cfg.chn_nse_grp{iGR},tme_rmv) > thr_rmv);
                        end
                        
                        fprintf('\n\nRemoving %i Channels\n\n',numel(unique(row)))
                        
                        rmv_chn{iGR} = [rmv_chn{iGR} ; dat.label(cfg.chn_nse_grp{iGR}(unique(row)))];
                        cfg.chn_nse_grp{iGR}(unique(row)) = [];
                        col_plt(unique(row),:) = [];
                        
                        close(ddd)
                        
                    else
                        
                        close(ddd)
                        
                        sat = 1;
                        
                        if thr_rmv > 5
                            
                            cfg.chn_nse_grp{rmv_hld_ind} = cellfun(@find,cellfun(@(x) strcmpi(ele_idn(:,1),x),rmv_chn{iGR},'UniformOutput',0));
                            cfg.pre_fix{rmv_hld_ind} = [cfg.pre_fix{iGR} '_2'];
                            
                            dat.cfg.rmv_nse.spl_chn{rmv_hld_ind} = cfg.chn_nse_grp{rmv_hld_ind};
                            dat.cfg.rmv_nse.spl_nme{rmv_hld_ind} = rmv_hld_nme;
                            dat.cfg.rmv_nse.pre_fix{rmv_hld_ind} = cfg.pre_fix{rmv_hld_ind};
                            
                        elseif tme_rmv > 5
                            
                            row = 1:numel(dat.label(cfg.chn_nse_grp{iGR}));
                            
                            rmv_chn{iGR} = [rmv_chn{iGR} ; dat.label(cfg.chn_nse_grp{iGR})];
                            cfg.chn_nse_grp{iGR}(unique(row)) = [];
                            col_plt(unique(row),:) = [];
                            
                        end
                    end
                    ttt_men = squeeze(mean(ttt_dat(cfg.chn_nse_grp{iGR},:,:),1))';
                end
                
            else
                rmv_chn{iGR} = {};
            end
            
            dat.cfg.cmn_nse_upd.(cfg.pre_fix{iGR}) = rmv_chn{iGR};
            if ~isempty(rmv_chn{iGR}) ; sve_num = rmv_chn{iGR}; else sve_num = {}; end
            save([cfg.out_dir '/' cfg.sbj_nme '/' cfg.pre_fix{iGR} '/'  'rmv_num.mat'],'sve_num')
        end
        
        if ~isempty(cfg.chn_nse_grp{iGR})
            if ~all(size(dat.trial{cfg.rnd_trl(iRN)}(cfg.chn_nse_grp{iGR}(1),:))==size(ttt_men(1,:))); ttt_men = ttt_men'; end
            
            for iCH = 1:numel(dat.label(cfg.chn_nse_grp{iGR}))
                for iTS = 1:numel(dat.trial)
                    
                    dat.trial{iTS}(cfg.chn_nse_grp{iGR}(iCH),:) = dat.trial{iTS}(cfg.chn_nse_grp{iGR}(iCH),:) - ttt_men(iTS,:);
                    
                end
            end
            
            ttt_dat = cat(3,dat.trial{:});
            ttt_men = squeeze(mean(ttt_dat(cfg.chn_nse_grp{iGR},:,:),1))';
            
            for iRN = 1:numel(cfg.rnd_tme)
                set(0,'currentfigure',h(iRN)); subplot(3,1,2); hold on;
                plot(dat.time{1},ttt_men(cfg.rnd_trl(iRN),:),'k');
                xlim([dat.time{1}(cfg.rnd_tme(iRN)) dat.time{1}(cfg.rnd_tme(iRN)+round(dat.fsample)-1)]);
                print(gcf,[cfg.out_dir '/' cfg.sbj_nme '/' cfg.pre_fix{iGR} '/' cfg.pre_fix{iGR} '_' num2str(iRN) '.png'],'-dpng');
            end
            
            for iRN = 1:numel(cfg.rnd_tme)
                set(0,'currentfigure',h(iRN)); subplot(3,1,3); hold on;
                for iP = 1:numel(dat.label(cfg.chn_nse_grp{iGR})); plot(dat.time{1},dat.trial{cfg.rnd_trl(iRN)}(cfg.chn_nse_grp{iGR}(iP),:),'Color',col_plt(iP,:)); end
                plot(dat.time{1},ttt_men(cfg.rnd_trl(iRN),:),'k','LineWidth',1.25);
                xlim([dat.time{1}(cfg.rnd_tme(iRN)) dat.time{1}(cfg.rnd_tme(iRN)+round(dat.fsample)-1)]);
                print(gcf,[cfg.out_dir '/' cfg.sbj_nme '/' cfg.pre_fix{iGR} '/' cfg.pre_fix{iGR} '_' num2str(iRN) '.png'],'-dpng');
            end
            
            close all
            
        end
        
    end
    
    iGR = iGR + 1;
    
end

end