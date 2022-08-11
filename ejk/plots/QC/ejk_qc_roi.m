function ejk_qc_roi(cfg)

if ~isfield(cfg,'zlm_plt'); cfg.zlm_plt = 3; end

%% Setup
ejk_chk_dir(cfg.out_dir)

col_hld = distinguishable_colors(size(cfg.dta,1));

%% Z-Score Data
% Z-Score By ROI
roi_zsc = nan(size(cfg.dta));
for iC = 1:size(cfg.dta,2)
    avg_dta = nanmean(cfg.dta(:,iC));
    dev_dta = nanstd(cfg.dta(:,iC));
    top_dta(1,iC) = (dev_dta*3)+avg_dta;
    bot_dta(1,iC) = (dev_dta*-3)+avg_dta;
    for iR = 1:size(cfg.dta,1)
        roi_zsc(iR,iC) = (cfg.dta(iR,iC) - avg_dta) / dev_dta;
    end
end

prb_ind_roi_zsc = [cfg.sbj_nme num2cell(sum(roi_zsc>2 | roi_zsc<-2,2)) num2cell(sum(roi_zsc>3 | roi_zsc<-3,2)) num2cell(sum(roi_zsc>4 | roi_zsc<-4,2))];

% Z-Score By Subj
sbj_zsc = nan(size(cfg.dta));
for iR = 1:size(cfg.dta,1)
    avg_dta = nanmean(cfg.dta(iR,:));
    dev_dta = nanstd(cfg.dta(iR,:));    
    for iC = 1:size(cfg.dta,2)
        sbj_zsc(iR,iC) = (cfg.dta(iR,iC) - avg_dta) / dev_dta;
    end
end

% Variance by subject
sbj_var = nan(size(cfg.dta,1),1);
for iR = 1:size(cfg.dta,1)
    sbj_var(iR,1) = nanstd(cfg.dta(iR,:));
end

% IQR by ROI
roi_iqr = nan(size(cfg.dta));
for iC = 1:size(cfg.dta,2)
    
    iqr_hld = iqr(cfg.dta(:,iC));
    qan_hld = quantile(cfg.dta(:,iC),[0.25 0.75]);
    
    for iR = 1:size(cfg.dta,1)
        if cfg.dta(iR,iC) > (qan_hld(2)+(iqr_hld*3)) || cfg.dta(iR,iC) < (qan_hld(1)-(iqr_hld*3))
            roi_iqr(iR,iC) = 2;
        elseif  cfg.dta(iR,iC) > (qan_hld(2)+(iqr_hld*1.5)) || cfg.dta(iR,iC) < (qan_hld(1)-(iqr_hld*1.5))
            roi_iqr(iR,iC) = 1;
        end
    end
end
prb_ind_roi_iqr = [ cfg.sbj_nme num2cell(sum(roi_iqr==1,2)) num2cell(sum(roi_iqr==2,2))];

[~, srt_ind] = sort(sbj_var);
sbj_var      = [prb_ind_roi_zsc(srt_ind,:) num2cell(sbj_var(srt_ind)) prb_ind_roi_iqr(srt_ind,:)];

%% Make Plots
% Raw %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot just raw
ylm = [ repmat(min(cfg.dta(:)), 1, size(cfg.dta,2)) ; repmat(max(cfg.dta(:)), 1, size(cfg.dta,2)) ];

figure('Position',[0 0 1080 1080],'Visible','off')
spider_plot( cfg.dta, 'Color', col_hld,'AxesLimits', ylm, ...
             'AxesColor', rgb('light grey')-0.1, ...
             'AxesDisplay', 'one', 'AxesInterval', 2, 'AxesLabels', cfg.dta_lbl, ...
             'MarkerSize', 12, 'MarkerSize', 1, 'LineWidth', 0.2 )
print(gcf,[cfg.out_dir '/' cfg.out_pre_fix '_' 'raw.png'],'-dpng')
close all

% Setup
sbj_inc =  find(sum(roi_zsc > cfg.zlm_plt | roi_zsc < -cfg.zlm_plt, 2));

raw_dta = [ top_dta ; bot_dta ; cfg.dta(sbj_inc,:)];
zsc_col = [ rgb('grey')-.15; rgb('grey')-.15; col_hld(sbj_inc,:)];
zne_wdt = [ 3, 3, repmat(1.5,1,numel(sbj_inc))];

% Plot 
ylm = [ repmat(min(raw_dta(:)), 1, size(raw_dta,2)) ; repmat(max(raw_dta(:)), 1, size(raw_dta,2)) ];

figure('Position',[0 0 1080 1080],'Visible','off')
spider_plot( raw_dta, 'Color', zsc_col,'AxesLimits', ylm, ...
             'AxesColor', rgb('light grey')-0.1, ...
             'AxesDisplay', 'one', 'AxesInterval', 2, 'AxesLabels', cfg.dta_lbl, ...
             'MarkerSize', 12, 'MarkerSize', 1, 'LineWidth', zne_wdt )
         
legend( cellfun(@(x) mmil_spec_char(x,{'_'}),[{'top'};{'bottom'};cfg.sbj_nme(sbj_inc)],'uni',0), 'Location', 'north','NumColumns',6, 'FontSize',12)
print(gcf,[cfg.out_dir '/' cfg.out_pre_fix '_' 'legend.png'],'-dpng')
close all

% Z by subject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ylm = [ repmat(min(sbj_zsc(:)), 1, size(sbj_zsc,2)) ; repmat(max(sbj_zsc(:)), 1, size(sbj_zsc,2)) ];

figure('Position',[0 0 1080 1080],'Visible','off')
spider_plot( sbj_zsc, 'Color', col_hld, 'AxesLimits', ylm, ...
             'AxesColor', rgb('light grey')-0.1, ...
             'AxesDisplay', 'one', 'AxesInterval', 2, 'AxesLabels', cfg.dta_lbl, ...
             'MarkerSize', 12, 'MarkerSize', 1, 'LineWidth', 0.2 )
print(gcf,[cfg.out_dir '/' cfg.out_pre_fix '_' 'zscore_by_subject.png'],'-dpng')
close all

% Z by region %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
ylm = [ repmat(min(roi_zsc(:)), 1, size(roi_zsc,2)) ; repmat(max(roi_zsc(:)), 1, size(roi_zsc,2)) ];

figure('Position',[0 0 1080 1080],'Visible','off')
spider_plot( roi_zsc, 'Color', col_hld, ...
             'AxesLimits', ylm, 'AxesColor', rgb('light grey')-0.1, ...
             'AxesDisplay', 'one', 'AxesInterval', 2, 'AxesLabels', cfg.dta_lbl, ...
             'MarkerSize', 1, 'LineWidth', 0.2 )
print(gcf,[cfg.out_dir '/' cfg.out_pre_fix '_' 'zscore_by_region.png'],'-dpng')

% Setup
sbj_inc =  find(sum(roi_zsc > cfg.zlm_plt | roi_zsc < -cfg.zlm_plt, 2));

zsc_dta = [ repmat(3,1,size(roi_zsc,2)) ; repmat(-3,1,size(roi_zsc,2)) ; roi_zsc(sbj_inc,:)];
zsc_col = [ rgb('dark grey')-.05; rgb('dark grey')-.05; col_hld(sbj_inc,:)];
zne_wdt = [ 3, 3, repmat(1.5,1,numel(sbj_inc))];

% Plot
ylm = [ repmat(min(zsc_dta(:)), 1, size(zsc_dta,2)) ; repmat(max(zsc_dta(:)), 1, size(zsc_dta,2)) ];

figure('Position',[0 0 1080 1080],'Visible','off')
spider_plot( zsc_dta, 'Color', zsc_col, ...
             'AxesLimits', ylm, 'AxesColor', rgb('light grey')-0.1, ...
             'AxesDisplay', 'one', 'AxesInterval', 2, 'AxesLabels', cfg.dta_lbl, ...
             'MarkerSize', 1, 'LineWidth', zne_wdt )
legend( cellfun(@(x) mmil_spec_char(x,{'_'}),[{'top'};{'bottom'};cfg.sbj_nme(sbj_inc)],'uni',0), 'Location', 'northoutside','NumColumns',6, 'FontSize',12)
print(gcf,[cfg.out_dir '/' cfg.out_pre_fix '_' 'zscore_by_region_legend.png'],'-dpng')
close all

%% Save Data
cell2csv( [cfg.out_dir '/' cfg.out_pre_fix '_' 'zscore_by_subject.csv']       , [ 'sbj_nme' cfg.dta_lbl ; cfg.sbj_nme num2cell(sbj_zsc)]);
cell2csv( [cfg.out_dir '/' cfg.out_pre_fix '_' 'zscore_by_region.csv']        , [ 'sbj_nme' cfg.dta_lbl ; cfg.sbj_nme num2cell(roi_zsc) ]);
cell2csv( [cfg.out_dir '/' cfg.out_pre_fix '_' 'zscore_by_region_zscores.csv'], prb_ind_roi_zsc);
cell2csv( [cfg.out_dir '/' cfg.out_pre_fix '_' 'subject_output.csv'] , sbj_var);

end