function [vox_dta , tot_lbl] = ejk_extract_fmri_roi_voxels(cfg)

%%
sbj_dir_lst = dir(sprintf('%s/FSURF*',cfg.fmr_fsr_dir));
n = regexp({sbj_dir_lst.name},['FSURF_' cfg.sbj_fsr_dir '.+_1$'],'match'); n = [n{:}];

if ~isempty(n)
    
    %% Get Labels
    if strcmpi(cfg.prc_nme,'') || strcmpi(cfg.prc_nme,'.a2009s')
        [~,fsr_lhs_lbl,~] = fs_read_annotation([cfg.fmr_fsr_dir '/' 'fsaverage' '/' 'label' '/' 'lh.aparc' cfg.prc_nme '.annot']);
        fsr_lhs_lbl = strcat('ctx-lh-',unique(fsr_lhs_lbl));
        [~,fsr_rhs_lbl,~] = fs_read_annotation([cfg.fmr_fsr_dir '/' 'fsaverage' '/' 'label' '/' 'rh.aparc' cfg.prc_nme '.annot']);
        fsr_rhs_lbl = strcat('ctx-rh-',unique(fsr_rhs_lbl));
        tot_lbl = [ fsr_lhs_lbl ; fsr_rhs_lbl ];
    end
    
    %% Create Stats
    
    
    %%
    if strcmpi(cfg.prc_nme,'') || strcmpi(cfg.prc_nme,'.a2009s')
        
        roi_nme_hld = mmil_readtext([cfg.fmr_fsr_dir '/' n{1} '/' 'label' '/' 'aparc.annot' cfg.prc_nme '.ctab'],'  ');
        roi_nme_hld = [roi_nme_hld(1:10,2:3) ; roi_nme_hld(11:end,1:2)];
        
        % DATA LOCATION
        if     exist([cfg.fmr_dir '/' cfg.sbj_fmr_dir '/' 'orig' '/' cfg.fmr_nme '_' 'cs20+orig.BRIK'])==2
            org_loc = 'orig';
            dta_loc = 'orig';
        elseif exist([cfg.fmr_dir '/' cfg.sbj_fmr_dir '/' 'orig.BLOCK' '/' cfg.fmr_nme '_' 'cs20+orig.BRIK'])==2
            org_loc = 'orig';
            dta_loc = 'orig.BLOCK';
        else
            dta_loc = [];
        end
        
        % RUN AFNI PART
        if ~isempty(dta_loc) && (exist('dta_loc') || exist([cfg.fmr_dir '/' cfg.sbj_fmr_dir '/' dta_loc]) ~= 7)
                        
            for iR = 1:size(roi_nme_hld,1)
                
                cmd = '';
                cmd     = sprintf('%scd %s/%s/%s\n',cmd,cfg.fmr_dir,cfg.sbj_fmr_dir,org_loc);
                cmd     = [cmd '3dROIstats -nomeanout ' cfg.fmr_typ ' -mask ' ...
                    'aparc' cfg.prc_nme '+aseg_rank_Alnd_Exp::ctx_lh_' roi_nme_hld{iR,2} ' ' ...
                    '../' dta_loc '/' cfg.fmr_nme '_' 'cs20+orig.BRIK' ' > ' ...
                    cfg.prj_dir '/' 'DATA' '/' cfg.sbj_nme '/' 'BOLD' '/' 'ROI' '/' 'aparc' '_' mmil_spec_char(cfg.prc_nme,{'.' '-' '&' ' '}) '_' mmil_spec_char(cfg.fmr_typ,{'.' '-' '&' ' '}) '_lh_' roi_nme_hld{iR,2} '_' cfg.sbj_nme '.1D'];
                out = unix(cmd);
                
                try lhs_hld = mmil_readtext([cfg.prj_dir '/' 'DATA' '/' cfg.sbj_nme '/' 'BOLD' '/' 'ROI' '/' 'aparc' '_' mmil_spec_char(cfg.prc_nme,{'.' '-' '&' ' '}) '_' mmil_spec_char(cfg.fmr_typ,{'.' '-' '&' ' '}) '_lh_' roi_nme_hld{iR,2} '_' cfg.sbj_nme '.1D'],'\t');
                    lhs_out(iR,:) = {['lh.' roi_nme_hld{iR,2}] lhs_hld{2,3}}; catch; end
                delete([cfg.prj_dir '/' 'DATA' '/' cfg.sbj_nme '/' 'BOLD' '/' 'ROI' '/' 'aparc' '_' mmil_spec_char(cfg.prc_nme,{'.' '-' '&' ' '}) '_' mmil_spec_char(cfg.fmr_typ,{'.' '-' '&' ' '}) '_lh_' roi_nme_hld{iR,2} '_' cfg.sbj_nme '.1D']);
                
                cmd = '';
                cmd     = sprintf('%scd %s/%s/%s\n',cmd,cfg.fmr_dir,cfg.sbj_fmr_dir,org_loc);
                cmd     = [cmd '3dROIstats -nomeanout ' cfg.fmr_typ ' -mask ' ...
                    'aparc' cfg.prc_nme '+aseg_rank_Alnd_Exp::ctx_rh_' roi_nme_hld{iR,2} ' ' ...
                    '../' dta_loc '/' cfg.fmr_nme '_' 'cs20+orig.BRIK' ' > ' ...
                    cfg.prj_dir '/' 'DATA' '/' cfg.sbj_nme '/' 'BOLD' '/' 'ROI' '/' 'aparc' '_' mmil_spec_char(cfg.prc_nme,{'.' '-' '&' ' '}) '_' mmil_spec_char(cfg.fmr_typ,{'.' '-' '&' ' '}) '_rh_' roi_nme_hld{iR,2} '_' cfg.sbj_nme '.1D'];
                unix(cmd);
                
                try rhs_hld = mmil_readtext([cfg.prj_dir '/' 'DATA' '/' cfg.sbj_nme '/' 'BOLD' '/' 'ROI' '/' 'aparc' '_' mmil_spec_char(cfg.prc_nme,{'.' '-' '&' ' '}) '_' mmil_spec_char(cfg.fmr_typ,{'.' '-' '&' ' '}) '_rh_' roi_nme_hld{iR,2} '_' cfg.sbj_nme '.1D'],'\t');
                    rhs_out(iR,:) = {['rh.' roi_nme_hld{iR,2}] rhs_hld{2,3}}; catch; end
                delete([cfg.prj_dir '/' 'DATA' '/' cfg.sbj_nme '/' 'BOLD' '/' 'ROI' '/' 'aparc' '_' mmil_spec_char(cfg.prc_nme,{'.' '-' '&' ' '}) '_' mmil_spec_char(cfg.fmr_typ,{'.' '-' '&' ' '}) '_rh_' roi_nme_hld{iR,2} '_' cfg.sbj_nme '.1D']);
                
            end
            
            lhs_out(cellfun(@isempty,lhs_out(:,1)),:) = [];
            rhs_out(cellfun(@isempty,rhs_out(:,1)),:) = [];
            
            vox_dta = [lhs_out ;  rhs_out];
            
            fprintf('Finished Subject : %s\n\n',cfg.sbj_nme)
            
        else
            
            fprintf('No Data For Subject : %s\n\n',cfg.sbj_nme);
            vox_dta = [];

        end
        
    else
        
        fprintf('No Data For Subject : %s\n\n',cfg.sbj_nme)
        vox_dta = [];
        
    end
    
else
    
    fprintf('No Data For Subject : %s\n\n',cfg.sbj_nme)
    vox_dta = [];
    tot_lbl = [];
    
end

end