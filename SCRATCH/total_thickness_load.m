clear; clc;

prj_dir = '/home/ekaestne/PROJECTS/';

prc_nme = { '' '.a2009s' '.split' };

%%
sbj_nme = mmil_readtext( [ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_subject_name.csv' ] );
fsr_nme = mmil_readtext( [ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_freesurfer_recons.csv' ] );

ovr_thk_dta = cell( size(sbj_nme,1) , 3);

%%
for iS = 1:size(sbj_nme,1)
    
    ovr_thk_dta{iS,1} = sbj_nme{iS,1};
    
    % Get Subject Data
    cfg.sbj_nme     = fsr_nme{iS,1};
    cfg.ovr_dir     = fsr_nme{iS,2};
    cfg.sbj_fsr_dir = fsr_nme{iS,3};
    
    sbj_dir_lst = dir(sprintf('%s/FSURF_*',cfg.ovr_dir));
    sbj_dir_lst = regexp({sbj_dir_lst.name},['FSURF_' cfg.sbj_fsr_dir '.+_1$'],'match'); sbj_dir_lst = [sbj_dir_lst{:}];
    
    % Load 
    if ~isempty(sbj_dir_lst) && numel(sbj_dir_lst)==1
        
        try
            
            lhs_ovr_dta = mmil_readtext([ cfg.ovr_dir '/' sbj_dir_lst{1} '/' 'stats' '/' 'lh.aparc.stats' ]);
            rhs_ovr_dta = mmil_readtext([ cfg.ovr_dir '/' sbj_dir_lst{1} '/' 'stats' '/' 'rh.aparc.stats' ]);
            
            ovr_thk_dta{iS,2} = lhs_ovr_dta{strcmpi(lhs_ovr_dta(:,3),'Mean Thickness'),4};
            ovr_thk_dta{iS,3} = rhs_ovr_dta{strcmpi(rhs_ovr_dta(:,3),'Mean Thickness'),4};
            
        catch
            ovr_thk_dta{iS,2} = nan;
            ovr_thk_dta{iS,3} = nan;
        end
    else
        ovr_thk_dta{iS,2} = nan;
        ovr_thk_dta{iS,3} = nan;
    end
    
end

%%
cell2csv([ '/home/ekaestne/PROJECTS/DATA/ROIHOLD/' '/' 'total_thickness.csv'], [ {'SubjID' 'LHS' 'RHS'} ; ovr_thk_dta])