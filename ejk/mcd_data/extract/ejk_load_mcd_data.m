
% cfg = [];
% cfg.prj_dir = prj_dir;
% cfg.dta_typ = 'fib_tfa'; % 'fib_tmd' 'fmr_lng_des_num_vox' 'fmr_lng_dst_num_vox' 'fmr_lng_ali_num_vox' 
% cfg.sbj_nme = grp_fle(:,1);
% 
% 'fib_tfa'             'fib_tmd'
% 'wmp_wfa_des'         'wmp_wmd_des'
% 'vol_dta'
% 'fmr_lng_des_num_vox' 'fmr_lng_dst_num_vox' 'fmr_lng_ali_num_vox'
% 

function dta = ejk_load_mcd_data(cfg)

%%
if isfield( cfg , 'fle_nme' )
    
    lod_dta_hld = mmil_readtext([ cfg.fle_nme ]);
    
else
    
    switch cfg.dta_typ
        case 'fib_tfa'
            lod_dta_hld = mmil_readtext([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'Fibers_aparc_FA.csv' ]);
        case 'fib_tmd'
            lod_dta_hld = mmil_readtext([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'Fibers_aparc_MD.csv' ]);
        case 'wmp_wfa_des'
            lod_dta_hld = mmil_readtext([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'WMParc_aparc__FA.csv' ]);
        case 'wmp_wmd_des'
            lod_dta_hld = mmil_readtext([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'WMParc_aparc__MD.csv' ]);
        case 'vol_dta'
            lod_dta_hld = mmil_readtext([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'Volumes_aparc.csv' ]);
        case 'gry_thk_des'
            lod_dta_hld = mmil_readtext([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'MRI_thickness_aparc_.csv' ]);
        case 'gry_thk_dst'
            lod_dta_hld = mmil_readtext([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'MRI_thickness_aparc_xa2009s.csv' ]);
        case 'fmr_lng_des_num_vox'
            
        case 'fmr_lng_dst_num_vox'
            lod_dta_hld = mmil_readtext([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'fMRI_aparc_xa2009s_N_FF_xnzvoxels.csv' ]);
        case 'fmr_lng_ali_num_vox'
            lod_dta_hld = mmil_readtext([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'fMRI_aparc_xAlicia_N_FF_xnzvoxels.csv' ]);
        case 'dti_cnn'
            lod_dta_hld = load([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'DTI_Connectome_norm.mat' ]);
    end
    
end

%%
if iscell(lod_dta_hld)
    sbj_nme_hld = lod_dta_hld(2:end,1);
    col_nme_hld = lod_dta_hld(1,2:end);
    dta_hld     = cell2mat(lod_dta_hld(2:end,2:end));
elseif isstruct(lod_dta_hld)
    sbj_nme_hld = lod_dta_hld.sve_sbj_nme;

    switch cfg.dta_typ
        case 'dti_cnn'
            col_nme_hld = lod_dta_hld.row_lbl;
            row_nme_hld = lod_dta_hld.col_lbl;
            dta_hld     = lod_dta_hld.cnn_dta;
    end
end

%%
dta.sbj_nme = sbj_nme_hld;

if iscell(lod_dta_hld)  
    for iCL = 1:numel(col_nme_hld)
        dta.(['x' mmil_spec_char(col_nme_hld{iCL},{'-'})]) = dta_hld(:,iCL);
    end
elseif isstruct(lod_dta_hld)
    switch cfg.dta_typ
        case 'dti_cnn'
            dta.dti_cnn_dta = dta_hld;
            dta.row_lbl = row_nme_hld;
            dta.col_lbl = col_nme_hld;
            for iR = 1:numel(row_nme_hld)
                for iC = 1:numel(col_nme_hld)
                    dta.cll_lbl{iR,iC} = [row_nme_hld{iR} '==' col_nme_hld{iC}];
                end
            end
    end
end

%%
if isfield(cfg,'sbj_nme')
    
    if strcmpi(cfg.sbj_nme{1},'SubjID'); cfg.sbj_nme(1) = []; end
    sbj_ind = ismember( sbj_nme_hld , cfg.sbj_nme );
    for iSI = 1:sum(sbj_ind)
            sbj_hld(iSI,1) = find(strcmpi(sbj_nme_hld , cfg.sbj_nme{iSI}),1);
    end
    sbj_ind = sbj_hld;
    
    dta.sbj_nme = dta.sbj_nme(sbj_ind,1);
    if iscell(lod_dta_hld)
        for iCL = 1:numel(col_nme_hld)
            dta.(['x' mmil_spec_char(col_nme_hld{iCL},{'-'})]) = dta.(['x' mmil_spec_char(col_nme_hld{iCL},{'-'})])(sbj_ind,1);
        end
    elseif isstruct(lod_dta_hld)
        switch cfg.dta_typ
            case 'dti_cnn'
                dta.dti_cnn_dta = dta.dti_cnn_dta(sbj_ind,:,:);
        end
    end
end

end