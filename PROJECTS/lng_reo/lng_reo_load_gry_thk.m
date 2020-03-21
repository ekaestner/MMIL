clear; clc;

%% Directories
prj_dir     = '/home/ekaestne/PROJECTS/';
prc_dir     = '/home/ekaestne/PROJECTS/DATA/'; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_mri ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf
mri_fsr_dir = '/home/mmilmcd/data/FSRECONS/';

prj_nme     = 'LanguageReorganization';
sbj_grp_nme = 'Lng_Reo_subject_list.csv';

%% Types
lbl_nme = { '' '.a2009s' };

%% Put together
sbj_grp = mmil_readtext([ prj_dir '/' 'SUBJECTS' '/' 'tmp' '/' prj_nme '/' sbj_grp_nme ]);
sbj_fsr = mmil_readtext([ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_freesurfer_names.csv' ]);
sbj_bld = mmil_readtext([ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_bold_names.csv' ]);

%%
tic;
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.ovr_dir = mri_fsr_dir; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf

for iPR = 1:numel(lbl_nme)
    
    fcfg.prc_nme = lbl_nme{iPR};
    
    for iS = 1:size(sbj_grp,1)
        
        fprintf( [ 'Subject #' num2str(iS) ' of ' num2str(size(sbj_grp,1)) ' : '  sbj_grp{iS,1} '\n' ] );
        
        sbj_fsr_ind = find( strcmpi( sbj_fsr(:,2) , sbj_grp{iS,1} ) );
        
        fcfg.sbj_fsr_dir = sbj_fsr{sbj_fsr_ind,3};
        fcfg.sbj_nme     = sbj_fsr{sbj_fsr_ind,2};
        
        [ gry_thk_dta_hld , tot_lbl ]= ejk_extract_grey_thickness(fcfg);
        if ~isempty(gry_thk_dta_hld); gry_thk_dta(iS,:) = gry_thk_dta_hld; tot_lbl_hld = tot_lbl; else gry_thk_dta(iS,:) = nan(1,size(gry_thk_dta,2)); end
        
    end
    
    sve_sbj_nme = sbj_fsr(:,2);
    
    ejk_chk_dir([ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'GreyThick' '/'])
    cell2csv(   [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'GreyThick' '/' 'GreyThick' '_' 'aparc' '_' mmil_spec_char(fcfg.prc_nme,{'.'}) '.csv'] , [ ['SubjId' ; sve_sbj_nme(1:size(gry_thk_dta,1))] [tot_lbl_hld ; num2cell(gry_thk_dta)] ] );
    
    clear gry_thk_dta
    
end
toc