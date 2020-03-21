function [pst_cog_scr , out_cat] = ejk_post_cognitive(cfg,cog_dta,sbj_scn)

if ~isfield(cfg,'neg_oly'); cfg.neg_oly = 0; end

%% Get fieldnames
thr_hld = [-1.5 1.5];

cog_fld_nme = fieldnames(cog_dta); cog_fld_nme(strcmpi(cog_fld_nme,'sbj_nme')) = [];
scn_fld_nme = fieldnames(sbj_scn); scn_fld_nme(strcmpi(scn_fld_nme,'sbj_nme')) = [];

raw_pst_cog_fld_nme = cog_fld_nme(string_find(cog_fld_nme,'raw.*_pst'));
raw_pre_cog_fld_nme = cellfun(@(x) x(1:end-4),raw_pst_cog_fld_nme,'uni',0);
nor_pst_cog_fld_nme = cog_fld_nme(string_find(cog_fld_nme,'nor.*_pst'));
nor_pre_cog_fld_nme = cellfun(@(x) x(1:end-4),nor_pst_cog_fld_nme,'uni',0);

con_ind = string_find(cog_dta.sbj_nme,'fc');

%% Post Score Calculations
pst_cog_scr.sbj_nme = cog_dta.sbj_nme;

for iFN = 1:numel(raw_pre_cog_fld_nme)
    switch raw_pre_cog_fld_nme{iFN}
        case 'log_mem_raw_scr_one'
            pst_cog_scr.log_mem_scr_one = ( ( cog_dta.log_mem_nor_scr_one_pst - cog_dta.log_mem_nor_scr_one ) - (12.1-10.2) ) / 2.14;
        case 'log_mem_raw_scr_two'
            pst_cog_scr.log_mem_scr_two = ( ( cog_dta.log_mem_nor_scr_two_pst - cog_dta.log_mem_nor_scr_two ) - (12.5-10.2) ) / 2.07;
        case 'cvl_lfr_raw_scr'
            pst_cog_scr.cvl_lfr_scr = ( ( cog_dta.cvl_lfr_raw_scr_pst - cog_dta.cvl_lfr_raw_scr ) - (11.74-10.19) ) / 1.86;
        case 'vp1_raw_scr'
            pst_cog_scr.vp1_scr = ( ( cog_dta.vp1_nor_scr_pst - cog_dta.vp1_nor_scr ) - (11.8-10.4) ) / 1.94;
        case 'vp2_raw_scr'
            pst_cog_scr.vp2_scr = ( ( cog_dta.vp2_nor_scr_pst - cog_dta.vp2_nor_scr ) - (11.1-10.5) ) / 1.88;
        case 'bnt_raw_scr'
            pst_cog_scr.bnt_scr = ( ( cog_dta.bnt_nor_scr_pst - cog_dta.bnt_nor_scr ) - (50.06 - 48.78) ) / 2.67;
        case 'ant_mem_raw_scr'
            pst_cog_scr.ant_mem_scr = ( ( cog_dta.ant_mem_raw_scr_pst - cog_dta.ant_mem_raw_scr ) - (49.42-49.08) ) / 1.09;
        case 'cat_flu_raw_scr'
            pst_cog_scr.cat_flu_scr = ( ( cog_dta.cat_flu_nor_scr_pst - cog_dta.cat_flu_nor_scr ) - (10.3-9.83) ) / 2.15;            
        case 'cvl_tot_raw_scr'
            pst_cog_scr.cvl_tot_scr = ( ( cog_dta.cvl_tot_raw_scr_pst - cog_dta.cvl_tot_raw_scr ) - (56.05-47.92) ) / 7.65;
        case 'ltr_tot_raw_scr'
            pst_cog_scr.ltr_tot_scr = ( ( cog_dta.ltr_tot_nor_scr_pst - cog_dta.ltr_tot_nor_scr ) - (10.1-9.62) ) / 2.11;
        case 'swt_cor_raw_scr'
            pst_cog_scr.swt_cor_scr = ( ( cog_dta.swt_cor_nor_scr_pst - cog_dta.swt_cor_nor_scr ) - (9.87-9.86) ) / 3.42;
        case 'swt_acc_raw_scr'
            pst_cog_scr.swt_acc_scr = ( ( cog_dta.swt_acc_nor_scr_pst - cog_dta.swt_acc_nor_scr ) - (10.54-10.41) ) / 3.75;
    end
end

%% Categorization
pst_cog_cat.sbj_nme = cog_dta.sbj_nme;

pst_fld_nme = fieldnames(pst_cog_scr); pst_fld_nme(1) = [];

if cfg.neg_oly
    for iFN = 1:numel(pst_fld_nme)
        
        pst_cog_cat.(pst_fld_nme{iFN}) = nan(numel(pst_cog_scr.(pst_fld_nme{iFN})),1);
        
        imp_ind = [pst_cog_scr.(pst_fld_nme{iFN}) < thr_hld(1)];
        nor_ind = [pst_cog_scr.(pst_fld_nme{iFN}) > thr_hld(1)];
        
        pst_cog_cat.(pst_fld_nme{iFN})(imp_ind) = 1;
        pst_cog_cat.(pst_fld_nme{iFN})(nor_ind) = 2;
        
    end
else
    for iFN = 1:numel(pst_fld_nme)
        
        pst_cog_cat.(pst_fld_nme{iFN}) = nan(numel(pst_cog_scr.(pst_fld_nme{iFN})),1);
        
        imp_ind = [pst_cog_scr.(pst_fld_nme{iFN}) < thr_hld(1)];
        nor_ind = [pst_cog_scr.(pst_fld_nme{iFN}) > thr_hld(1) & pst_cog_scr.(pst_fld_nme{iFN}) < thr_hld(2)];
        inc_ind = [pst_cog_scr.(pst_fld_nme{iFN}) > thr_hld(2)];
        
        pst_cog_cat.(pst_fld_nme{iFN})(imp_ind) = 1;
        pst_cog_cat.(pst_fld_nme{iFN})(nor_ind) = 2;
        pst_cog_cat.(pst_fld_nme{iFN})(inc_ind) = 3;
        
    end
end

%% Plots
% if ~exist([cfg.out_dir '/' 'cognitive']); mkdir([cfg.out_dir '/' 'cognitive']); end
% 
% % Get number of subjects with scan
% for iSC = 1:2:numel(scn_fld_nme)
%     sbj_has_hld{iSC} = find(sbj_scn.(scn_fld_nme{iSC}));
%     for iT = 1:numel(pst_fld_nme)
%         sbj_has_pst(iT) = numel(intersect(find(~isnan(pst_cog_cat.(pst_fld_nme{iT}))),sbj_has_hld{iSC}));
%         sbj_pst_imp(iT) = numel(intersect(find(pst_cog_cat.(pst_fld_nme{iT})==1),sbj_has_hld{iSC}));
%         sbj_pst_nrm(iT) = numel(intersect(find(pst_cog_cat.(pst_fld_nme{iT})==2),sbj_has_hld{iSC}));
%         sbj_pst_sup(iT) = numel(intersect(find(pst_cog_cat.(pst_fld_nme{iT})==3),sbj_has_hld{iSC}));
%     end
%     tbl = [ [{''}  pst_fld_nme'] ; [{'Total post-op patients'} num2cell(sbj_has_pst)] ; [{'Impaired'} num2cell(sbj_pst_imp)] ; [{'No Change'} num2cell(sbj_pst_nrm)] ; [{'Improvement'} num2cell(sbj_pst_sup)] ];
%     
%     cell2csv([cfg.out_dir '/' 'cognitive' '/' scn_fld_nme{iSC} '.csv'],tbl)
%     
% end
% 
% %
% for iT = 1:numel(pst_fld_nme)
%     sbj_has_pst(iT) = numel(find(~isnan(pst_cog_cat.(pst_fld_nme{iT}))));
%     sbj_pst_imp(iT) = numel(find(pst_cog_cat.(pst_fld_nme{iT})==1));
%     sbj_pst_nrm(iT) = numel(find(pst_cog_cat.(pst_fld_nme{iT})==2));
%     sbj_pst_sup(iT) = numel(find(pst_cog_cat.(pst_fld_nme{iT})==3));
% end
% tbl = [ [{''}  pst_fld_nme'] ; [{'Total post-op patients'} num2cell(sbj_has_pst)] ; [{'Impaired'} num2cell(sbj_pst_imp)] ; [{'No Change'} num2cell(sbj_pst_nrm)] ; [{'Improvement'} num2cell(sbj_pst_sup)] ];
% 
% cell2csv([cfg.out_dir '/' 'cognitive' '/' scn_fld_nme{iSC} '.csv'],tbl)
% 
% % Plot of available subjects
% sub_plt_dim = ceil(sqrt(numel(pst_fld_nme)));
% 
% figure('units','normalized','outerposition',[0 0 1 1],'Visible','off')
% for iT = 1:numel(pst_fld_nme)
%     
%     con_sbj{iT}     = string_find(cog_dta.sbj_nme,'fc');
%     epd_sbj_pre{iT} = string_find(cog_dta.sbj_nme,'epd');
%     
%     epd_sbj_pst_imp{iT} = find(pst_cog_cat.(pst_fld_nme{iT})==1);
%     epd_sbj_pst_nrm{iT} = find(pst_cog_cat.(pst_fld_nme{iT})==2);
%     epd_sbj_pst_sup{iT} = find(pst_cog_cat.(pst_fld_nme{iT})==3);
%     epd_sbj_pst_non{iT} = find(isnan(pst_cog_cat.(pst_fld_nme{iT})));
%            
%     xvl_con = (0.6 - 0.4) .* rand(numel(con_sbj{iT}),1) + 0.4;
%     xvl_pre = (0.9 - 0.7) .* rand(size(cog_dta.sbj_nme,1),1) + 0.7;
%     
%     subplot(sub_plt_dim,sub_plt_dim,iT); hold on;
%     scatter(xvl_con,cog_dta.(raw_pre_cog_fld_nme{iT})(con_sbj{iT},1),'k','filled')
%     
%     scatter(xvl_pre(epd_sbj_pst_non{iT}),cog_dta.(raw_pre_cog_fld_nme{iT})(epd_sbj_pst_non{iT},1),'k')
%     scatter(xvl_pre(epd_sbj_pst_nrm{iT}),cog_dta.(raw_pre_cog_fld_nme{iT})(epd_sbj_pst_nrm{iT},1),'b','filled')
%     scatter(xvl_pre(epd_sbj_pst_imp{iT}),cog_dta.(raw_pre_cog_fld_nme{iT})(epd_sbj_pst_imp{iT},1),'r','filled')
%     scatter(xvl_pre(epd_sbj_pst_sup{iT}),cog_dta.(raw_pre_cog_fld_nme{iT})(epd_sbj_pst_sup{iT},1),'g','filled')
%         
%    for iP = 1:numel(epd_sbj_pst_nrm{iT})
%         plot([1 1.5],[cog_dta.(raw_pre_cog_fld_nme{iT})(epd_sbj_pst_nrm{iT}(iP),1) cog_dta.(raw_pst_cog_fld_nme{iT})(epd_sbj_pst_nrm{iT}(iP),1)],'b','LineWidth',2); end
%      for iP = 1:numel(epd_sbj_pst_imp{iT})
%         plot([1 1.5],[cog_dta.(raw_pre_cog_fld_nme{iT})(epd_sbj_pst_imp{iT}(iP),1) cog_dta.(raw_pst_cog_fld_nme{iT})(epd_sbj_pst_imp{iT}(iP),1)],'r','LineWidth',2); end
%     for iP = 1:numel(epd_sbj_pst_sup{iT})
%         plot([1 1.5],[cog_dta.(raw_pre_cog_fld_nme{iT})(epd_sbj_pst_sup{iT}(iP),1) cog_dta.(raw_pst_cog_fld_nme{iT})(epd_sbj_pst_sup{iT}(iP),1)],'g','LineWidth',2); end
%     title(upper(mmil_spec_char(pst_fld_nme{iT}(1:end-4),{'_'})))
%     xlim([0.1 1.8])
%     
% end
% tightfig();
% 
% print([cfg.out_dir '/' 'cognitive' '/' 'overall_tasks' '.png'],'-dpng')

%% Put together
% if ~exist([cfg.out_dir '/' 'cognitive' '/' 'subjects' '/']); mkdir([cfg.out_dir '/' 'cognitive' '/' 'subjects' '/']); end

lrn_nme = { 'log_mem_scr_one' 'cvl_tot_scr' 'vp1_scr'};
mem_nme = { 'log_mem_scr_two' 'cvl_lfr_scr' 'vp2_scr' };
lng_nme = { 'bnt_scr' 'ant_mem_scr' 'cat_flu_scr' };
frn_nme = { 'ltr_tot_scr' 'swt_cor_scr' 'swt_acc_scr' };

ovr_nme = [lrn_nme {''} mem_nme {''} lng_nme {''} frn_nme ];

for iFN = 1:numel(ovr_nme)
    if ~strcmpi(ovr_nme{iFN},'')
        for iS = 1:size(pst_cog_scr.(ovr_nme{iFN}),1)
            dta_hld(iS,iFN) = pst_cog_scr.(ovr_nme{iFN})(iS);
        end
    else
        for iS = 1:size(pst_cog_scr.(ovr_nme{iFN-1}),1)
            dta_hld(iS,iFN) = nan;
        end
    end
    tot_num(:,iFN) = ~isnan(dta_hld(:,iFN));
end
tot_num = sum(max(tot_num'));

dta_hld(dta_hld>3) = 3;
dta_hld(dta_hld<-3) = -3;

dst_col = distinguishable_colors(tot_num);

%% Plot
% col_num = 1;
% for iS = 1:size(dta_hld,1)
%     if any(~isnan(dta_hld(iS,:)))
%         
%         figure('Visible','off'); hold on;
%         plot(dta_hld(iS,:),'LineWidth',3,'Color',dst_col(col_num,:),'marker','o','markeredgecolor',rgb('black'));
%         col_num = col_num + 1;
%         
%         xlim([0.5 size(dta_hld,2)+0.5])
%         
%         for iL = 1:size(dta_hld,2)
%             line([iL-.5 iL-.5],[-4 4],'LineWidth',0.5,'LineStyle','--','Color',rgb('light grey')-[0.15 0.15 0.15])
%         end
%         
%         hline(-1.5,'r')
%         hline(1.5,'g')
%         ylim([-4 4])
%         
%         xticks(1:size(dta_hld,2))
%         set(gca,'XTickLabel',cellfun(@(x) mmil_spec_char(x,{'_'}),ovr_nme,'uni',0))
%         xtickangle(60)
%         
%         tightfig();
%         
%         print([cfg.out_dir '/' 'cognitive' '/' 'subjects' '/' cog_dta.sbj_nme{iS} '_tasks' '.png'],'-dpng')
%         close all
%         
%     end
% end

%% Output
cat_nme = { 'Impaired' 'NoChange' 'Improved' };

out_scr = nan(size(dta_hld,1),sum(~strcmpi(ovr_nme,'')));

out_cat = repmat({''},size(dta_hld,1),sum(~strcmpi(ovr_nme,'')));

% Cognitive Put Together
col_num = 1;
for iC = 1:size(ovr_nme,2)
    if ~strcmpi(ovr_nme{iC},'')
        for iS = 1:size(dta_hld,1)
            if ~isnan(pst_cog_scr.(ovr_nme{iC})(iS))
                
                out_scr(iS,col_num) = pst_cog_scr.(ovr_nme{iC})(iS);
                
                out_cat{iS,col_num} = cat_nme{pst_cog_cat.(ovr_nme{iC})(iS)};
                                
            end
        end
        col_num = col_num+1;
    end
end

out_scr = [ovr_nme(~strcmpi(ovr_nme,'')) ; num2cell(out_scr)];
out_cat = [ovr_nme(~strcmpi(ovr_nme,'')) ; out_cat];

out_scr = [[{''} ; cog_dta.sbj_nme] out_scr];
out_cat = [[{''} ; cog_dta.sbj_nme] out_cat];

% Save
% cell2csv([cfg.out_dir '/' 'cognitive' '/' 'post_cognitive_scores.csv'],out_scr);
% cell2csv([cfg.out_dir '/' 'cognitive' '/' 'post_cognitive_categories.csv'],out_cat);

%% Put together Tables
% if ~exist([cfg.out_dir '/' 'cognitive' '/' 'tables' '/']); mkdir([cfg.out_dir '/' 'cognitive' '/' 'tables' '/']); end
% 
% cat_nme = { 'Impaired' 'NoChange' 'Improved' };
% 
% pst_cog_cat_nme = fieldnames(pst_cog_cat);
% 
% out_lab = {'Subject Name' 'Category' 'RCI' '' 'Side Onset' 'MTS' 'Age of Onset' '# AEDs' 'Seizure Frequency' '' 'Age' 'Sex' 'Handedness' 'Education' };
% 
% for iF = 1:numel(pst_cog_cat_nme)
%     
%     out_tab = cell(sum(~isnan(pst_cog_cat.(pst_cog_cat_nme{iF}))),14);
%     
%     psb_cat = unique(pst_cog_cat.(pst_cog_cat_nme{iF})(~isnan(pst_cog_cat.(pst_cog_cat_nme{iF}))));
%     for iPS = 1:numel(psb_cat)
%         
%         cat_loc = find(pst_cog_cat.(pst_cog_cat_nme{iF})==psb_cat(iPS));
%         
%         out_tab(find(cellfun(@isempty,out_tab(:,1)),1):find(cellfun(@isempty,out_tab(:,1)),1)+numel(cat_loc)-1,:) = ...
%             [ cog_dta.sbj_nme(cat_loc) ...
%               repmat(cat_nme(psb_cat(iPS)),numel(cat_loc),1) ...
%               num2cell(pst_cog_scr.(pst_cog_cat_nme{iF})(cat_loc)) ...
%               ...
%               repmat({''},numel(cat_loc),1) ...
%               sbj_sze.sbj_sde_ons(cat_loc) ...
%               sbj_sze.sbj_mts(cat_loc) ...
%               num2cell(sbj_sze.sbj_age_ons(cat_loc)) ...
%               num2cell(sbj_sze.sbj_aed_num(cat_loc)) ...
%               num2cell(sbj_sze.sbj_sze_frq(cat_loc)) ...
%               ...
%               repmat({''},numel(cat_loc),1) ...
%               sbj_dem.sbj_age(cat_loc) ...
%               sbj_dem.sbj_sex(cat_loc) ...
%               sbj_dem.sbj_hnd(cat_loc) ...
%               num2cell(sbj_dem.sbj_edu(cat_loc)) ...
%             ];
%         
%     end
%     
%     % Sort table
%    for iPS = 1:numel(psb_cat)
%         
%         ovr_srt_ind = find(strcmpi(out_tab(:,2),cat_nme{psb_cat(iPS)}));
%                 
%         % MTS - L
%         lft_ons_lft_mts_ind = find(strcmpi(out_tab(ovr_srt_ind,5),'L') & strcmpi(out_tab(ovr_srt_ind,6),'L'));
%         lft_ons_rgh_mts_ind = find(strcmpi(out_tab(ovr_srt_ind,5),'L') & strcmpi(out_tab(ovr_srt_ind,6),'R'));
%         lft_ons_non_mts_ind = find(strcmpi(out_tab(ovr_srt_ind,5),'L') & (strcmpi(out_tab(ovr_srt_ind,6),'N/A') | strcmpi(out_tab(ovr_srt_ind,6),'')));
%         
%         % MTS - R
%         rgh_ons_lft_mts_ind = find(strcmpi(out_tab(ovr_srt_ind,5),'R') & strcmpi(out_tab(ovr_srt_ind,6),'L'));
%         rgh_ons_rgh_mts_ind = find(strcmpi(out_tab(ovr_srt_ind,5),'R') & strcmpi(out_tab(ovr_srt_ind,6),'R'));
%         rgh_ons_non_mts_ind = find(strcmpi(out_tab(ovr_srt_ind,5),'R') & (strcmpi(out_tab(ovr_srt_ind,6),'N/A') | strcmpi(out_tab(ovr_srt_ind,6),'')));
%         
%         % MTS - ???
%         non_ons_lft_mts_ind = find((strcmpi(out_tab(ovr_srt_ind,5),'N/A') | strcmpi(out_tab(ovr_srt_ind,5),'')) & strcmpi(out_tab(ovr_srt_ind,6),'L'));
%         non_ons_rgh_mts_ind = find((strcmpi(out_tab(ovr_srt_ind,5),'N/A') | strcmpi(out_tab(ovr_srt_ind,5),'')) & strcmpi(out_tab(ovr_srt_ind,6),'R'));
%         non_ons_non_mts_ind = find((strcmpi(out_tab(ovr_srt_ind,5),'N/A') | strcmpi(out_tab(ovr_srt_ind,5),'')) & (strcmpi(out_tab(ovr_srt_ind,6),'N/A') | strcmpi(out_tab(ovr_srt_ind,6),'')));
%         
%         if numel(ovr_srt_ind) ~= numel([ lft_ons_lft_mts_ind ; lft_ons_rgh_mts_ind ; lft_ons_non_mts_ind ; rgh_ons_rgh_mts_ind ; rgh_ons_lft_mts_ind ; rgh_ons_non_mts_ind ; non_ons_rgh_mts_ind ; non_ons_lft_mts_ind ; non_ons_non_mts_ind ]); error(':p'); end
%         
%         out_tab(ovr_srt_ind,:) = out_tab(ovr_srt_ind([ lft_ons_lft_mts_ind ; lft_ons_rgh_mts_ind ; lft_ons_non_mts_ind ; rgh_ons_rgh_mts_ind ; rgh_ons_lft_mts_ind ; rgh_ons_non_mts_ind ; non_ons_rgh_mts_ind ; non_ons_lft_mts_ind ; non_ons_non_mts_ind ]),:);
%         
%     end
%     
%     % Out    
%     out_tab = [out_lab ; out_tab];
%     
%     cell2csv([cfg.out_dir '/' 'cognitive' '/' 'tables' '/' pst_cog_cat_nme{iF} '.csv'],out_tab);
%     
% end
% 
% 
% end

% %%
% pst_fld_nme
% raw_pst_cog_fld_nme
% raw_pre_cog_fld_nme
% nor_pst_cog_fld_nme
% nor_pre_cog_fld_nme
% 
% iF = 12;
% [num2cell([ pst_cog_cat.(pst_fld_nme{iF}) pst_cog_scr.(pst_fld_nme{iF}) cog_dta.(raw_pre_cog_fld_nme{iF}) cog_dta.(raw_pst_cog_fld_nme{iF}) cog_dta.(nor_pre_cog_fld_nme{iF-1}) cog_dta.(nor_pst_cog_fld_nme{iF-1})]) cog_dta.sbj_nme]






