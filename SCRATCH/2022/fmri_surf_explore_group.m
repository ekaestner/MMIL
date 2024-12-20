clear; clc;

%% Inputs
% Locations
sbj_dir = '/space/mcdonald-syn01/1/data/MCD_MRI/UWFaces_bidsified/derivatives/afni/';

bld_dir = '/space/mcdonald-syn01/1/data/MCD_MRI/UWFaces_bidsified/derivatives/afni/GroupAnalysis/';

sma_dir = '/home/mmilmcd2/.afni/data/suma_MNI152_2009';

out_dir = '/home/ekaestne/PROJECTS/OUTPUT/fce_nme_emy/';

% To use
fnc_mri_stt = {'All_Tstat' };
fnc_mri_nme = {'all' };
fnc_mri_srf = {'pial' 'inflated'  };

sma_fle_nme = 'MNI152_2009';

fnc_mri_typ = {'-nzvoxels'};

% Parameters
fmr_rng_num = [-7 7];
cls_sze = 20;
sig_thr = 2.59;

% Stats
bld_fle_nme = { 'New_GLT' 'New-Repeated_GLT' };
bld_fnc_nme = { 'eventTent' 'block' };

%% Constants
fmr_col_map = {'neon blue' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

hms_nme = { 'lhs' 'rhs' };

%% Subject find
sbj_nme = {'Ttestpp'};
grp_cmp = { 'full_group' 'HC' 'LTLE' 'RTLE' 'LTLEvsRTLE' 'TLE' 'TLEvsHC' };

%% Make the Stat Mask
for iGC = 1:numel(grp_cmp)
    for iB = 1:numel(bld_fle_nme)
        for iF = 1:numel(bld_fnc_nme)
            for iST = 1:numel(fnc_mri_nme)
                for iS = 1:numel(sbj_nme)
                    
                    ejk_chk_dir([ out_dir '/' sbj_nme{iS} '/']);
                    ejk_chk_dir([ out_dir '/' sbj_nme{iS} '/' 'stats' '/']);
                    ejk_chk_dir([ out_dir '/' sbj_nme{iS} '/' 'stats' '/' grp_cmp{iGC} '/' ]);
                    
                    fcfg = [];
                    
                    fcfg.sbj_nme = sbj_nme{iS};
                    
                    fcfg.bld_dir     = [ bld_dir '/' sbj_nme{iS} '/' ];
                    fcfg.bld_fle_nme = [ grp_cmp{iGC} '_desc-3dttestpp_tmap_' bld_fle_nme{iB} '#0_Coef' '_' bld_fnc_nme{iF} ];
                    
                    fcfg.fnc_mri_stt = fnc_mri_stt{iST};
                    fcfg.fnc_mri_nme = fnc_mri_nme{iST};
                    fcfg.mdl_fnc_nme = bld_fnc_nme{iF};
                    
                    fcfg.cls_sze = cls_sze;
                    fcfg.sig_thr = sig_thr;
                    
                    fcfg.out_dir = [ out_dir '/' sbj_nme{iS} '/' 'stats' '/' grp_cmp{iGC} '/'];
                    fcfg.pre_fix = [ fcfg.sbj_nme '_' bld_fle_nme{iB} '_' fcfg.fnc_mri_nme '_' 'Tstat' '_' fcfg.mdl_fnc_nme ];
                    
                    ejk_fmri_stat_create_nii(fcfg)
                    
                end
            end
        end
    end
end

%% Pull The Surface Out
for iGC = 1:numel(grp_cmp)
    for iB = 1:numel(bld_fle_nme)
        for iF = 1:numel(bld_fnc_nme)
            for iST = 1:numel(fnc_mri_nme)
                for iS = 1:numel(sbj_nme)
                    
                    ejk_chk_dir([ out_dir '/' sbj_nme{iS} '/' 'surf' '/' ]);
                    ejk_chk_dir([ out_dir '/' sbj_nme{iS} '/' 'surf' '/' grp_cmp{iGC} '/']);
                    
                    iSU = find(strcmpi(fnc_mri_srf,'pial'));
                    
                    fcfg = [];
                    
                    fcfg.sbj_nme = sbj_nme{iS};
                    
                    fcfg.bld_dir = [ out_dir '/' sbj_nme{iS} '/' 'stats' '/' grp_cmp{iGC} '/' ];
                    fcfg.sma_dir = sma_dir;
                    
                    fcfg.bld_fle_nme = [ grp_cmp{iGC} '_desc-3dttestpp_tmap_' bld_fle_nme{iB} '#0_Coef' '_' bld_fnc_nme{iF} ];
                    fcfg.sma_fle_nme = sma_fle_nme;
                    
                    fcfg.fnc_mri_stt = fnc_mri_stt{iST};
                    fcfg.fnc_mri_nme = fnc_mri_nme{iST};
                    fcfg.fnc_mri_srf = fnc_mri_srf{iSU};
                    fcfg.mdl_fnc_nme = bld_fnc_nme{iF};
                    
                    fcfg.grd_par = [ fcfg.bld_dir '/' fcfg.sbj_nme '_' bld_fle_nme{iB} '_' fcfg.fnc_mri_nme '_' 'Tstat' '_' fcfg.mdl_fnc_nme ];;
                    
                    fcfg.cls_sze = cls_sze;
                    fcfg.sig_thr = sig_thr;
                    
                    fcfg.out_dir = [ out_dir '/' sbj_nme{iS} '/' 'surf' '/' grp_cmp{iGC} '/' ];
                    fcfg.pre_fix = [ fcfg.out_dir '/' fcfg.sbj_nme '_' bld_fle_nme{iB} '_' fcfg.fnc_mri_nme '_' 'Tstat' '_' fcfg.fnc_mri_srf '_' fcfg.mdl_fnc_nme ];
                    
                    ejk_fmri_extract_surf_nii(fcfg)
                    
                end
            end
        end
    end
end

%% Make fMRI surface Plot
for iGC = 1:numel(grp_cmp)
    for iB = 1:numel(bld_fle_nme)
        for iF = 1:numel(bld_fnc_nme)
            for iST = 1:numel(fnc_mri_nme)
                for iS = 1:numel(sbj_nme)
                    
                    fprintf([ 'plotting ' grp_cmp{iGC} ' ' '(' num2str(iGC) ' of ' num2str(numel(grp_cmp)) ')' ' | ' bld_fle_nme{iB} ' ' '(' num2str(iB) ' of ' num2str(numel(bld_fle_nme)) ')' ' | ' sbj_nme{iS} ' (' num2str(iS) ' of ' num2str(numel(sbj_nme)) ')' '\n']);
                    
                    ejk_chk_dir([ out_dir '/' sbj_nme{iS} '/' 'plot' '/' grp_cmp{iGC} '/']);

                    for iH = 1:numel(hms_nme)                        
                        % Load
                        fmr_val{iH} = afni_niml_read([ out_dir '/' sbj_nme{iS} '/' 'surf' '/' grp_cmp{iGC} '/' sbj_nme{iS} '_' bld_fle_nme{iB} '_' fnc_mri_nme{iST} '_' 'Tstat' '_' fnc_mri_srf{1} '_' bld_fnc_nme{iF} '_' hms_nme{iH} '.niml.dset']);
                        fmr_val{iH} = fmr_val{iH}{1}.nodes{1}.data; 
                    end
 
                    %
                    for iSU = 1:numel(fnc_mri_srf)
                        
                        % Both directions
                        pcfg = [];
                        
                        pcfg.fsr_dir = sma_dir;
                        
                        pcfg.srf_typ = fnc_mri_srf{iSU};
                        
                        pcfg.out_dir     = [ out_dir '/' sbj_nme{iS} '/' 'plot' '/' grp_cmp{iGC} '/'];
                        pcfg.out_pre_fix = [ sbj_nme{iS} '_' bld_fle_nme{iB} '_' fnc_mri_nme{iST} '_' fnc_mri_srf{iSU} '_' bld_fnc_nme{iF} '_' 'both'];
                        
                        pcfg.plt_dta = { fmr_val{1}(:,end)' fmr_val{2}(:,end)' };
                        
                        pcfg.fmr_col_map = {'neon blue' 'blue' 'greyish blue' 'grey' 'reddish grey' 'red' 'reddish orange' };
                        
                        pcfg.low_rng_num = [ -sig_thr sig_thr ];
                        pcfg.hgh_rng_num = fmr_rng_num ;
                        
                        mmil_anat_surf_plot(pcfg)
                        
                        % Positive only
                        pcfg.plt_dta{1}(pcfg.plt_dta{1}<0) = 0;
                        pcfg.plt_dta{2}(pcfg.plt_dta{2}<0) = 0;
                        
                        pcfg.out_pre_fix = [ sbj_nme{iS} '_' bld_fle_nme{iB} '_' fnc_mri_nme{iST} '_' fnc_mri_srf{iSU} '_' bld_fnc_nme{iF} '_' 'positive'];
                        
                        mmil_anat_surf_plot(pcfg)
                        
                    end
                end
            end
        end
    end
end
