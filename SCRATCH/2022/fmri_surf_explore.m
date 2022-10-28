clear; clc;

%% Inputs
% Locations
sbj_dir = '/space/mcdonald-syn01/1/data/MCD_MRI/UWFaces_bidsified/derivatives/afni/';
sbj_idn = 'sub-';

bld_dir = [ sbj_dir '/SBJNME/ses-1/func/'];

sma_dir = '/home/mmilmcd2/.afni/data/suma_MNI152_2009';

out_dir = '/home/ekaestne/PROJECTS/OUTPUT/fce_nme_emy/';

% To use
fnc_mri_stt = {'New_GLT#0'         'New-Repeated_GLT#0' };
fnc_mri_nme = {'Novel_vs_Baseline' 'Novel_vs_Repeated' };
fnc_mri_srf = {'pial' 'inflated'  };

sma_fle_nme = 'MNI152_2009';

fnc_mri_typ = {'-nzvoxels'};

% Parameters
fmr_rng_num = [-7 7];
cls_sze = 20;
sig_thr = 2.59;

% Stats
bld_fle_nme = 'SBJNME_ses-1_task-faces_space-MNI152NLin2009cAsym_desc-stats';
bld_fnc_nme = { 'eventTent' 'block' };

%% Constants
fmr_col_map = {'neon blue' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};

hms_nme = { 'lhs' 'rhs' };

%% Subject find
sbj_nme = dir(sbj_dir);
    sbj_nme = {sbj_nme([sbj_nme.isdir]).name};
    sbj_nme = sbj_nme(string_find(sbj_nme,sbj_idn));
    
%% Make the Stat Mask
for iF = 1:numel(bld_fnc_nme)
    for iST = 1:numel(fnc_mri_nme)
        for iS = 1:numel(sbj_nme)
            
            ejk_chk_dir([ out_dir '/' sbj_nme{iS} '/']);
            ejk_chk_dir([ out_dir '/' sbj_nme{iS} '/' 'stats' '/']);
            
            fcfg = [];
            
            fcfg.sbj_nme = sbj_nme{iS};
            
            fcfg.bld_dir     = [strrep(bld_dir,'SBJNME',sbj_nme{iS}) ];
            fcfg.bld_fle_nme = [strrep(bld_fle_nme,'SBJNME',sbj_nme{iS}) '_' bld_fnc_nme{iF}];
            
            fcfg.fnc_mri_stt = fnc_mri_stt{iST};
            fcfg.fnc_mri_nme = fnc_mri_nme{iST};
            fcfg.mdl_fnc_nme = bld_fnc_nme{iF}; 
            
            fcfg.cls_sze = cls_sze;
            fcfg.sig_thr = sig_thr;
            
            fcfg.out_dir = [ out_dir '/' sbj_nme{iS} '/' 'stats' '/'];
            
            ejk_fmri_stat_create_nii(fcfg)
            
        end
    end
end

%% Pull The Surface Out
for iF = 1:numel(bld_fnc_nme)
    for iST = 1:numel(fnc_mri_nme)
        for iS = 1:numel(sbj_nme)
            
            ejk_chk_dir([ out_dir '/' sbj_nme{iS} '/' 'surf' '/']);
            
            iSU = find(strcmpi(fnc_mri_srf,'pial'));
            
            fcfg = [];
            
            fcfg.sbj_nme = sbj_nme{iS};
            
            fcfg.bld_dir = [ out_dir '/' sbj_nme{iS} '/' 'stats' '/'];
            fcfg.sma_dir = sma_dir;
            
            fcfg.bld_fle_nme = strrep(bld_fle_nme,'SBJNME',sbj_nme{iS});
            fcfg.sma_fle_nme = sma_fle_nme;
            
            fcfg.fnc_mri_stt = fnc_mri_stt{iST};
            fcfg.fnc_mri_nme = fnc_mri_nme{iST};
            fcfg.fnc_mri_srf = fnc_mri_srf{iSU};
            fcfg.mdl_fnc_nme = bld_fnc_nme{iF};
            
            fcfg.cls_sze = cls_sze;
            fcfg.sig_thr = sig_thr;
            
            fcfg.out_dir = [ out_dir '/' sbj_nme{iS} '/' 'surf' '/'];
            
            ejk_fmri_extract_surf_nii(fcfg)
            
        end
    end
end

%% Make fMRI surface Plot
for iF = 1:numel(bld_fnc_nme)
    for iST = 1:numel(fnc_mri_nme)
        for iS = 1:numel(sbj_nme)
            
            ejk_chk_dir([ out_dir '/' sbj_nme{iS} '/' 'plot' '/']);
            fprintf([ 'plotting ' bld_fnc_nme{iF} ' ' '(' num2str(iF) ' of ' num2str(numel(bld_fnc_nme)) ')' ' | ' fnc_mri_nme{iST} ' ' '(' num2str(iST) ' of ' num2str(numel(fnc_mri_nme)) ')' ' | ' sbj_nme{iS} ' (' num2str(iS) ' of ' num2str(numel(sbj_nme)) ')' '\n']);
                    
            for iH = 1:numel(hms_nme)                
                % Load
                fmr_val{iH} = afni_niml_read([ out_dir '/' sbj_nme{iS} '/' 'surf' '/' sbj_nme{iS} '_' fnc_mri_nme{iST} '_' fnc_mri_srf{1} '_' 'Tstat' '_' bld_fnc_nme{iF} '_' hms_nme{iH} '.niml.dset']);
                fmr_val{iH} = fmr_val{iH}{1}.nodes{1}.data;                 
            end

            %
            for iSU = 1:numel(fnc_mri_srf)
                
                % Both directions
                pcfg = [];
                
                pcfg.fsr_dir = sma_dir;
                
                pcfg.srf_typ = fnc_mri_srf{iSU};
                
                pcfg.out_dir     = [ out_dir '/' sbj_nme{iS} '/' 'plot' '/'];
                pcfg.out_pre_fix = [ sbj_nme{iS} '_' fnc_mri_nme{iST} '_' fnc_mri_srf{iSU} '_' bld_fnc_nme{iF} '_' 'both'];
                
                pcfg.plt_dta = { fmr_val{1}(:,end)' fmr_val{2}(:,end)' };
                
                pcfg.fmr_col_map = {'neon blue' 'blue' 'greyish blue' 'grey' 'reddish grey' 'red' 'reddish orange' };
                
                pcfg.low_rng_num = [ -sig_thr sig_thr ];
                pcfg.hgh_rng_num = fmr_rng_num ;
                
                mmil_anat_surf_plot(pcfg)
                
                % Positive only
                pcfg.plt_dta{1}(pcfg.plt_dta{1}<0) = 0;
                pcfg.plt_dta{2}(pcfg.plt_dta{2}<0) = 0;
                
                pcfg.out_pre_fix = [ sbj_nme{iS} '_' fnc_mri_nme{iST} '_' fnc_mri_srf{iSU} '_' bld_fnc_nme{iF} '_' 'positive'];
                
                mmil_anat_surf_plot(pcfg)
                
            end
        end
    end
end
