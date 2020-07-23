% cfg = [];
%
% cfg.srf_dir = '/space/syn09/1/data/MMILDB/MCD_RSI/fsurf';
% cfg.prc_dir = '/space/syn09/1/data/MMILDB/MCD_RSI/proc_dti';
% cfg.sbj_nme = {'fc096' 'epd082'};
% cfg.sbj_dir = {'fc096_fmri_170728_20170728.171240_1' 'epd082_fmri2_160720_20160720.161406_1'};
%
% cfg.prc_nme = '';
%
% cfg.min_val = 1e-06;

function [ wmp_avg , col_nme ] = ejk_extract_wmparc(cfg)

if ~isfield(cfg,'min_val'); cfg.min_val = 1e-06; end

%%
sbj_dir_lst = dir(sprintf('%s/DTIPROC_*',cfg.prc_dir));
sbj_dir_lst = regexp({sbj_dir_lst.name},['DTIPROC_' cfg.sbj_fsr_dir '.+_1$'],'match'); sbj_dir_lst = [sbj_dir_lst{:}];

col_loc = mmil_readtext(['/home/ekaestne/PROJECTS/SCRIPTS/streamint' '/' 'wmparc' '.aparc' cfg.prc_nme '.' 'annot']);

if ~isempty(sbj_dir_lst)
    
    try
        
        vol_dta_nme = [cfg.prc_dir '/' sbj_dir_lst{1} '/' 'DTcalc' '/' 'DTI1_crev_corr_regT1_DT_' cfg.fib_typ '.mgz'];
        vol_dta_sze = ejk_read_vol_sze(vol_dta_nme,1);
        vol_dta = ejk_load_vol(vol_dta_nme);
        
        
        seg_dta_nme = [cfg.prc_dir '/' sbj_dir_lst{1} '/' 'DTanalysis' '/' 'wmparc_resDTI.mgz']; % - EJK - Fix to allow multiple label types
        seg_dta_sze = ejk_read_vol_sze(seg_dta_nme,1);
        seg_dta = ejk_load_vol(seg_dta_nme);
        
        for iFC = 1:size(col_loc,1)
            
            val_hld = vol_dta(find(seg_dta==col_loc{iFC,1})) * cfg.scl_fct;
            num_val = size(val_hld,1);
            
            ind_vld = find(min(abs(val_hld),[],2)>=cfg.min_val & ~isnan(sum(val_hld,2)));
            ind_invalid = setdiff([1:num_val],ind_vld);
            ind_nan = isnan(sum(val_hld,2));
            
            val_hld(ind_nan) = 0;
            
            if isempty(val_hld)
                wmp_avg(1,iFC) = nan;
                wmp_med(1,iFC) = nan;
                wmp_std(1,iFC) = nan;
            else
                wmp_avg(1,iFC) = nanmean(val_hld(ind_vld,:),1);
                wmp_med(1,iFC) = nanmedian(val_hld(ind_vld,:),1);
                if numel(ind_vld)>1
                    wmp_std(1,iFC) = nanstd(val_hld(ind_vld,:),1);
                end
            end
            
            col_nme(1,iFC) = col_loc(iFC,2);
            
        end
        
    catch
        for iFC = 1:size(col_loc,1)
            col_nme(1,iFC) = col_loc(iFC,2);
        end
        
        wmp_avg     = [];
        
    end
    
else
    for iFC = 1:size(col_loc,1)
        col_nme(1,iFC) = col_loc(iFC,2);
    end
    
    wmp_avg     = [];
    
end

end
