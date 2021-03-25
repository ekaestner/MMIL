% cfg = [];
%
% cfg.ovr_dir = '/space/syn09/1/data/MMILDB/MCD_RSI/proc_dti';
% cfg.sbj_nme = {'fc096' 'epd082'};
% cfg.sbj_fsr_dir = {'fc096_fmri_170728_20170728.171240_1' 'epd082_fmri2_160720_20160720.161406_1'};
%
% cfg.atl_dir = 'AtlasTrack';
% cfg.map_dir = 'fiber_maps';
%
% cfg.thr_prb = 0.08;
%
% cfg.fib_cde = [101 102 103 104 105 106 107 108 109 110 115 116 117 118 119 120 121 122 123 133 134 135 136 137 138 141 142 143 144 145 146 147 148 149 150 1014 1024 2000 2001 2002 2003 2004];
%
% cfg.min_val = 1e-06;

% extra
%
% atl_nme = [atl_loc '/' cfg.atl_nme];
% cfg.atl_nme = 'T1_Atlas/T1_atlas';

function [ fib_avg , col_nme ] = ejk_extract_fibers(cfg)

%%
if ~isfield(cfg,'min_val'); cfg.min_val = 1e-6; end

atl_loc = [getenv('MMPS_DIR') '/' 'parms' '/' 'DTI_Fiber' '/' ];
atl_leg = mmil_readtext([atl_loc '/' 'DTI_Fiber_Legend.csv']);

sbj_dir_lst = dir(sprintf('%s/DTIPROC*',cfg.ovr_dir));
sbj_dir_lst = regexp({sbj_dir_lst.name},['DTIPROC_' cfg.sbj_fsr_dir '.+_1$'],'match'); sbj_dir_lst = [sbj_dir_lst{:}];

if ~isempty(sbj_dir_lst)
    
    vol_dta_nme = dir([cfg.ovr_dir '/' sbj_dir_lst{1} '/' 'DTcalc' '/' '*.mgz']);
    vol_dta_nme_ind = regexp({vol_dta_nme.name},['crev_corr_regT1_DT_' cfg.fib_typ '.mgz']);
    vol_dta_nme = [vol_dta_nme(~cellfun(@isempty,vol_dta_nme_ind)).name];
    
    if ~isempty(vol_dta_nme)
        
        try
            
            vol_dta_sze = ejk_read_vol_sze([cfg.ovr_dir '/' sbj_dir_lst{1} '/' 'DTcalc' '/' vol_dta_nme],1);
            vol_dta = ejk_load_vol([cfg.ovr_dir '/' sbj_dir_lst{1} '/' 'DTcalc' '/' vol_dta_nme]);
            
            for iFC = 1:numel(cfg.fib_cde)
                
                fib_fle_nme = [cfg.ovr_dir '/' sbj_dir_lst{1} '/' 'AtlasTrack' '/' 'fiber_maps' '/' 'fiber' '_' num2str(cfg.fib_cde(iFC)) '_' 'prob' '_' 'countatlas' '_' 'xcg' '_' 'pthresh' num2str(cfg.thr_prb) '.mat' ];
                roi_vol_sze = ejk_read_vol_sze(fib_fle_nme);
                vol_roi = ejk_load_vol(fib_fle_nme);
                roi_ind = find(vol_roi>0);
                
                fin_val = vol_dta(roi_ind) * cfg.scl_fct;
                
                num_val = size(fin_val,1);
                
                ind_vld = find(abs(fin_val(:,1))>=cfg.min_val & ~isnan(fin_val(:,1)));
                ind_inv = setdiff(1:num_val,ind_vld);
                ind_nan = isnan(fin_val(:,1));
                
                fin_val(ind_nan,:) = 0;
                
                fin_wgh = vol_roi(roi_ind);
                fin_wgh(ind_inv) = 0;
                
                fib_avg(1,iFC) = mmil_wtd_mean(fin_val(ind_vld,:),fin_wgh(ind_vld,:),1);
                fib_med(1,iFC) = mmil_wtd_median(fin_val(ind_vld,:),fin_wgh(ind_vld,:),1);
                if numel(ind_vld)>1
                    fib_std(1,iFC) = mmil_wtd_std(fin_val(ind_vld,:),fin_wgh(ind_vld,:),1);
                end
                
                col_nme(1,iFC) = atl_leg(find(cell2mat(atl_leg(2:end,1))==cfg.fib_cde(iFC))+1,2);
                
            end
            
        catch
            
            col_nme = atl_leg(ismember(cell2mat(atl_leg(2:end,1)),cfg.fib_cde),2)';
            fib_avg = nan(1,size(col_nme,2));
            
        end
        
    else
        
        col_nme = atl_leg(ismember(cell2mat(atl_leg(2:end,1)),cfg.fib_cde),2)';
        fib_avg = nan(1,size(col_nme,2));
        
    end
    
else
    
    col_nme = atl_leg(ismember(cell2mat(atl_leg(2:end,1)),cfg.fib_cde),2)';
    fib_avg = nan(1,size(col_nme,2));
    
end
