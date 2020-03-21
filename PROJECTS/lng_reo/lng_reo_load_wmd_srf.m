clear; clc;

tic;

%% Directories
prj_dir     = '/home/ekaestne/PROJECTS/';
prc_dir     = '/home/ekaestne/PROJECTS/DATA/'; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf
dti_fsr_dir = '/home/mmilmcdRSI/data/proc_dti/';

prj_nme = 'LanguageReorganization';
sbj_grp_nme = 'Lng_Reo_subject_list.csv';

%% Types
lbl_nme = { '' '.a2009s' };

fib_typ = { 'FA' 'MD' }; % 'FA' 'MD'
scl_fct = [ 1    1000 ];
prc_nme = { '' '.a2009s' }; % '' '.a2009s'

%% Put together
sbj_grp = mmil_readtext([ prj_dir '/' 'SUBJECTS' '/' 'tmp' '/' prj_nme '/' sbj_grp_nme ]);
sbj_fsr = mmil_readtext([ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_freesurfer_names.csv' ]);
sbj_bld = mmil_readtext([ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_bold_names.csv' ]);

%%
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prc_dir = dti_fsr_dir; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf

for iPR = 1:numel(cfg.lbl_nme)
    
    fcfg.prc_nme = prc_nme{iPR};
    
    for iFT = 1:numel(fib_typ)
        
        fcfg.fib_typ = fib_typ{iFT};
        fcfg.scl_fct = scl_fct(iFT);
        
        for iS = 1:size(sbj_grp,1) % 6 10 42 57 58 59 80 82 87 89 90
            
            fprintf( [ 'Subject #' num2str(iS) ' of ' num2str(size(sbj_grp,1)) ' : '  sbj_grp{iS,1} '\n' ] );
            
            sbj_fsr_ind = find( strcmpi( sbj_fsr(:,2) , sbj_grp{iS,1} ) );
            
            fcfg.sbj_fsr_dir = sbj_fsr{sbj_fsr_ind,3};
            fcfg.sbj_nme     = sbj_fsr{sbj_fsr_ind,2};
            
            [ wmp_dta_hld , tot_lbl ]= ejk_extract_wmparc(fcfg);
            
            if ~isempty(wmp_dta_hld); wmp_dta(iS,:) = wmp_dta_hld; else wmp_dta(iS,:) = nan(1,size(tot_lbl,2)); end
            
        end
        
        sve_sbj_nme = sbj_fsr(:,2);
        
        ejk_chk_dir([ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'WMParc' '/'])
        cell2csv([ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'WMParc' '/' 'WMParc' '_' 'aparc' '_' mmil_spec_char(fcfg.prc_nme,{'.'}) '_' fib_typ{iFT} '.csv'],[ ['SubjId' ; sve_sbj_nme(1:size(wmp_dta,1))] [tot_lbl ; num2cell(wmp_dta)] ]);
        
        clear wmp_dta
        
    end
end
toc

