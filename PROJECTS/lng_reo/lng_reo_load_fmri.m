clear; clc;

tic;

%% Directories
prj_dir     = '/home/ekaestne/PROJECTS/';
prc_dir     = '/home/ekaestne/PROJECTS/DATA/'; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf
fmr_dir     = '/home/mmilmcd/data/MCD_BOLD/subjects';
fmr_fsr_dir = '/home/mmilmcdRSI/data/fsurf/';

prj_nme = 'LanguageReorganization';
sbj_grp_nme = 'Lng_Reo_subject_list.csv';

%% Types
lbl_nme = { '' '.a2009s' };

fmr_stt = 'N-FF_GLT#0_Tstat';
fmr_nme = 'N-FF';
fmr_typ = '-nzvoxels';

%% Put together
sbj_grp = mmil_readtext([ prj_dir '/' 'SUBJECTS' '/' 'tmp' '/' prj_nme '/' sbj_grp_nme ]);
sbj_fsr = mmil_readtext([ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_freesurfer_names.csv' ]);
sbj_bld = mmil_readtext([ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_bold_names.csv' ]);

%%
for iP = 2 % 1:numel(lbl_nme)
    
    prc_nme = lbl_nme{iP};
    
%     fmr_dta = nan( size(sbj_grp,1) , [] );

    for iS = 1:size(sbj_grp,1) % 18 78 79 80 81 82 128 129
        
        fprintf( [ 'Subject #' num2str(iS) ' of ' num2str(size(sbj_grp,1)) ' : '  sbj_grp{iS,1} '\n' ] );
        
        sbj_fsr_ind = find( strcmpi( sbj_fsr(:,2) , sbj_grp{iS,1} ) );
        sbj_bld_ind = find( strcmpi( sbj_bld(:,2) , sbj_grp{iS,1} ) );
        
        if ~isempty( sbj_bld{sbj_bld_ind,3} )
            
            fcfg = [];
               
            fcfg.prj_dir  = prj_dir;
            fcfg.fmr_dir = fmr_dir;
            fcfg.fmr_fsr_dir = fmr_fsr_dir;
            
            fcfg.prc_nme = prc_nme;
            
            fcfg.fmr_stt = fmr_stt;
            fcfg.fmr_nme = fmr_nme;
            fcfg.fmr_typ = fmr_typ;
            
            fcfg.sbj_fmr_dir = sbj_bld{ sbj_bld_ind , 3 };
            fcfg.sbj_fsr_dir = sbj_fsr{ sbj_fsr_ind , 3 };
            fcfg.sbj_nme     = sbj_grp{ iS , 1 };
            
            [ fmr_dta_hld , tot_lbl ] = ejk_extract_fmri_roi_voxels(fcfg);
            
            if ~isempty(fmr_dta_hld)
                
                for iRW = 1:size(fmr_dta_hld(:,1),1)
                    
                    nme_ind = regexp( fmr_dta_hld{iRW,1} , '\.' , 'split' );
                    nme_ind = intersect( string_find( tot_lbl , nme_ind{1} ) , string_find( tot_lbl , nme_ind{2} ) );
                    
                    if ~isempty(nme_ind)
                        fmr_dta(iS,nme_ind) = fmr_dta_hld{ iRW , 2 };
                    else
                        fmr_dta(iS,nme_ind) = 0;
                    end
                    
                end
                
                prc_nme_hld = tot_lbl;
                
            else
                
                fmr_dta(iS,:) = nan(1,size(fmr_dta,2));
                
            end
        else
            
            fmr_dta(iS,:) = nan(1,size(fmr_dta,2));
            
        end
    end
    
    sve_sbj_nme = sbj_fsr(:,2);
    
    ejk_chk_dir( [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'fMRI' '/' ] )
    save(     [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'fMRI' '/' 'fMRI' '_' 'aparc' '_' mmil_spec_char(prc_nme,{'.'}) '_' mmil_spec_char(fmr_nme,{'.' '-' ' '}) '_' mmil_spec_char(fmr_typ,{'.' '-' ' '}) '.mat'],'sve_sbj_nme','prc_nme','fmr_dta');
    cell2csv( [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'fMRI' '/' 'fMRI' '_' 'aparc' '_' mmil_spec_char(prc_nme,{'.'}) '_' mmil_spec_char(fmr_nme,{'.' '-' ' '}) '_' mmil_spec_char(fmr_typ,{'.' '-' ' '}) '.csv'],[ ['SubjId' ; sve_sbj_nme(1:size(fmr_dta,1))] [prc_nme_hld' ; num2cell(fmr_dta)] ])
    
    clear fmr_dta
    
end
toc

