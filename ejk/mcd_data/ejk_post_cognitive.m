function pst_cog_scr = ejk_post_cognitive(cfg,cog_dta)

if ~isfield(cfg,'neg_oly'); cfg.neg_oly = 0; end
if ~isfield(cfg,'thr_hld'); thr_hld = [-1.5 1.5]; end
if ~isfield(cfg,'rci');     cfg.rci = 1; end

%% Get fieldnames
cog_fld_nme = fieldnames(cog_dta); cog_fld_nme(strcmpi(cog_fld_nme,'sbj_nme')) = [];

raw_pst_cog_fld_nme = cog_fld_nme(string_find(cog_fld_nme,'raw.*_pst'));
raw_pre_cog_fld_nme = cellfun(@(x) x(1:end-4),raw_pst_cog_fld_nme,'uni',0);
nor_pst_cog_fld_nme = cog_fld_nme(string_find(cog_fld_nme,'nor.*_pst'));
nor_pre_cog_fld_nme = cellfun(@(x) x(1:end-4),nor_pst_cog_fld_nme,'uni',0);

pst_nme = [ raw_pst_cog_fld_nme ; nor_pst_cog_fld_nme];
pre_nme = [ raw_pre_cog_fld_nme ; nor_pre_cog_fld_nme];

pst_cog_scr.sbj_nme = cog_dta.sbj_nme;

%% Raw Change Scores
for iFN = 1:numel(pst_nme)
    pst_cog_scr.raw.(pst_nme{iFN}) = cog_dta.(pst_nme{iFN}) - cog_dta.(pre_nme{iFN});
end

%% RCI score calculation
for iFN = 1:numel(pre_nme)
    switch pre_nme{iFN}
        case 'log_mem_raw_scr_one'
            for iS = 1:numel(cog_dta.log_mem_nor_scr_one)
                pst_cog_scr.rci.log_mem_nor_scr_one_pst(iS,1) = ( ( cog_dta.log_mem_nor_scr_one_pst(iS) - cog_dta.log_mem_nor_scr_one(iS) ) - (12.1-10.2) ) / 2.14;
            end
        case 'log_mem_raw_scr_two'
            for iS = 1:numel(cog_dta.log_mem_nor_scr_two)
                if strcmpi(cog_dta.wms_ver_log_mem{iS},'III')
                    pst_cog_scr.rci.log_mem_nor_scr_two_pst(iS,1) = ( ( cog_dta.log_mem_nor_scr_two_pst(iS) - cog_dta.log_mem_nor_scr_two(iS) ) - (12.5-10.2) ) / 2.07;
                elseif strcmpi(cog_dta.wms_ver_log_mem{iS},'IV')
                    pst_cog_scr.rci.log_mem_nor_scr_two_pst(iS,1) = ( ( cog_dta.log_mem_nor_scr_two_pst(iS) - cog_dta.log_mem_nor_scr_two(iS) ) - (12.6-10.3) ) / 2.17;
                else
                    pst_cog_scr.rci.log_mem_nor_scr_two_pst(iS,1) = NaN;
                end
            end
        case 'cvl_lfr_raw_scr'
            for iS = 1:numel(cog_dta.cvl_lfr_raw_scr)
                pst_cog_scr.rci.cvl_lfr_raw_scr_pst(iS,1)     = ( ( cog_dta.cvl_lfr_raw_scr_pst(iS) - cog_dta.cvl_lfr_raw_scr(iS) ) - (11.74-10.19) ) / 1.86;
            end
        case 'vp1_raw_scr'
            for iS = 1:numel(cog_dta.vp1_raw_scr)
                pst_cog_scr.rci.vp1_nor_scr_pst(iS,1)         = ( ( cog_dta.vp1_nor_scr_pst(iS) - cog_dta.vp1_nor_scr(iS) ) - (11.8-10.4) ) / 1.94;
            end
        case 'vp2_raw_scr'
            for iS = 1:numel(cog_dta.vp2_raw_scr)
                if strcmpi(cog_dta.wms_ver_vpa{iS},'II')
                    pst_cog_scr.rci.vp2_nor_scr_pst(iS,1) = NaN;
                elseif strcmpi(cog_dta.wms_ver_vpa{iS},'III')
                    pst_cog_scr.rci.vp2_nor_scr_pst(iS,1)         = ( ( cog_dta.vp2_nor_scr_pst(iS) - cog_dta.vp2_nor_scr(iS) ) - (11.1-10.5) ) / 1.88;
                elseif strcmpi(cog_dta.wms_ver_vpa{iS},'IV')
                    pst_cog_scr.rci.vp2_nor_scr_pst(iS,1)         = ( ( cog_dta.vp2_nor_scr_pst(iS) - cog_dta.vp2_nor_scr(iS) ) - (10.8-9.8) ) / 1.98;
                else
                    pst_cog_scr.rci.vp2_nor_scr_pst(iS,1) = NaN;
                end
            end
        case 'bnt_raw_scr'
            for iS = 1:numel(cog_dta.bnt_raw_scr)
                pst_cog_scr.rci.bnt_raw_scr_pst(iS,1)         = ( ( cog_dta.bnt_raw_scr_pst(iS) - cog_dta.bnt_raw_scr(iS) ) - (50.06 - 48.78) ) / 2.67;
            end
        case 'ant_mem_raw_scr'
            for iS = 1:numel(cog_dta.ant_mem_raw_scr)
                pst_cog_scr.rci.ant_mem_raw_scr_pst(iS,1)     = ( ( cog_dta.ant_mem_raw_scr_pst(iS) - cog_dta.ant_mem_raw_scr(iS) ) - (49.42-49.08) ) / 1.09;
            end
        case 'cat_flu_raw_scr'
            for iS = 1:numel(cog_dta.cat_flu_raw_scr)
                pst_cog_scr.rci.cat_flu_nor_scr_pst(iS,1)     = ( ( cog_dta.cat_flu_nor_scr_pst(iS) - cog_dta.cat_flu_nor_scr(iS) ) - (10.3-9.83) ) / 2.15;
            end
        case 'cvl_tot_raw_scr'
            for iS = 1:numel(cog_dta.cvl_tot_raw_scr)
                pst_cog_scr.rci.cvl_tot_raw_scr_pst(iS,1)     = ( ( cog_dta.cvl_tot_raw_scr_pst(iS) - cog_dta.cvl_tot_raw_scr(iS) ) - (56.05-47.92) ) / 7.65;
            end
        case 'ltr_tot_raw_scr'
            for iS = 1:numel(cog_dta.ltr_tot_raw_scr)
                pst_cog_scr.rci.ltr_tot_nor_scr_pst(iS,1)     = ( ( cog_dta.ltr_tot_nor_scr_pst(iS) - cog_dta.ltr_tot_nor_scr(iS) ) - (10.1-9.62) ) / 2.11;
            end
        case 'swt_cor_raw_scr'
            for iS = 1:numel(cog_dta.swt_cor_raw_scr)
                pst_cog_scr.rci.swt_cor_nor_scr_pst(iS,1)     = ( ( cog_dta.swt_cor_nor_scr_pst(iS) - cog_dta.swt_cor_nor_scr(iS) ) - (9.87-9.86) ) / 3.42;
            end
        case 'swt_acc_raw_scr'
            for iS = 1:numel(cog_dta.swt_acc_raw_scr)
                pst_cog_scr.rci.swt_acc_nor_scr_pst(iS,1)     = ( ( cog_dta.swt_acc_nor_scr_pst(iS) - cog_dta.swt_acc_nor_scr(iS) ) - (10.54-10.41) ) / 3.75;
            end
    end
end

%% Percent Change Calculation
for iFN = 1:numel(pst_nme)
    pst_cog_scr.pct.(pst_nme{iFN}) = round(((cog_dta.(pst_nme{iFN}) - cog_dta.(pre_nme{iFN})) ./ cog_dta.(pre_nme{iFN}))*100);
end


