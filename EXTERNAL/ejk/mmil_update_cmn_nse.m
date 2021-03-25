function tmp = mmil_update_cmn_nse(cfg,dat)

%% Update cmn_nse
if ~isfield(cfg,'spl_chn')
    cur_cmn_nse = mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'cmn_nse');
    cur_cmn_nse_nme = mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'cmn_nme');
else
    cur_cmn_nse = cfg.spl_chn;  
    cur_cmn_nse_nme = cfg.spl_nme;
    
    rmv_ind = cellfun(@isempty,cur_cmn_nse{1});
    cur_cmn_nse{1}(rmv_ind) = [];
    cur_cmn_nse_nme{1}(rmv_ind) = [];
    cfg.pre_fix(rmv_ind) = [];
end

cur_cmn_nse = cur_cmn_nse{cfg.sub_fld_ind};
hld_org_nse = cur_cmn_nse;
cur_cmn_nse_nme = cur_cmn_nse_nme{cfg.sub_fld_ind};

new_cmn_nse_rmv = cell(1);
new_cmn_nse = cell(1);
for iCN = 1:numel(cur_cmn_nse)   
      
    load([cfg.out_dir '/' cfg.pre_fix{iCN} '/' 'rmv_num.mat'])
    [~,new_cmn_nse_rmv{iCN},~] = intersect(dat.label(cur_cmn_nse{iCN}),sve_num);
    new_cmn_nse{iCN} = cur_cmn_nse{iCN};
    new_cmn_nse{iCN}(new_cmn_nse_rmv{iCN}) = [];
        
end

new_str_hld = '{';
for iN = 1:numel(new_cmn_nse)
    new_str_hld = [new_str_hld num_2_str(new_cmn_nse{iN})]; 
    new_str_hld = [new_str_hld ','];
end
new_str_hld(end) = '}';

if ~isfield(cfg,'spl_chn')
    cur_cmn_nse = mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'cmn_nse','asstring');
    cur_cmn_nse = strsplit(cur_cmn_nse{1},';');
    cur_cmn_nse{cfg.sub_fld_ind} = new_str_hld;
    cur_cmn_nse = [strjoin(cur_cmn_nse,';') ';'];
    mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'cmn_nse',cur_cmn_nse);
else
   cur_cmn_nse = new_str_hld;
   mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'cmn_nse',[cur_cmn_nse ';']);
   mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'cmn_nme',['{' strjoin(cur_cmn_nse_nme,',') '};']);
end

mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'nse_val',num2str(2));

%% Save plots about what it would look like to subtract the missing channels (and save their form for easy insertion if desired)
cfg.rnd_ind = sort(round(randsample(dat.time{1},5)));
ttt_dat = cat(3,dat.trial{:});
for iP = 1:numel(hld_org_nse)
    
    oth_chn{iP} = setxor(hld_org_nse{iP},new_cmn_nse{iP});
        
    if ~isempty(oth_chn{iP})
        
        ttt_men = squeeze(mean(ttt_dat(oth_chn{iP},:,:),1))';
        
        % 1st Plot
        if ~isdir([cfg.out_dir '/' cfg.pre_fix{iP} '/' 'rmv_chn']); mkdir([cfg.out_dir '/' cfg.pre_fix{iP} '/' 'rmv_chn']); end
        col_plt = distinguishable_colors(numel(oth_chn{iP}));
        for iRN = 1:numel(cfg.rnd_ind)
            h(iRN) = figure('Visible','off'); subplot(2,1,1); hold on;
            for iPL = 1:numel(dat.label(oth_chn{iP})); plot(dat.time{1},dat.trial{1}(oth_chn{iP}(iPL),:),'Color',col_plt(iPL,:)); end
            plot(dat.time{1},ttt_men,'k','LineWidth',1.25);
            xlim([cfg.rnd_ind(iRN) cfg.rnd_ind(iRN)+2]);
        end
        
        % Subtraction
        for iCH = 1:numel(dat.label(oth_chn{iP}))
            for iTS = 1:numel(dat.trial)
                
                dat.trial{iTS}(oth_chn{iP}(iCH),:) = dat.trial{iTS}(oth_chn{iP}(iCH),:) - ttt_men(iTS,:);
                
            end
        end
        
        % 2nd Plot
        for iRN = 1:numel(cfg.rnd_ind)
            set(0,'currentfigure',h(iRN)); subplot(2,1,2); hold on;
            for iPL = 1:numel(dat.label(oth_chn{iP})); plot(dat.time{1},dat.trial{1}(oth_chn{iP}(iPL),:),'Color',col_plt(iPL,:)); end
            plot(dat.time{1},ttt_men,'k','LineWidth',1.25);
            xlim([cfg.rnd_ind(iRN) cfg.rnd_ind(iRN)+2]);
            print(gcf,[cfg.out_dir '/' cfg.pre_fix{iP} '/' 'rmv_chn' '/' cfg.pre_fix{iP} '_' num2str(iRN) '.png'],'-dpng');
        end
        
        % Save #'s
        tmp_chn_str = num_2_str(oth_chn{iP});
        cell2csv([cfg.out_dir '/' cfg.pre_fix{iP} '/' 'rmv_chn' '/' 'chn_str.csv'],{tmp_chn_str})

    end
end

tmp = dat;

end
















