clear; clc;

tic;

%% Directories
prj_dir     = '/home/ekaestne/PROJECTS/';
prc_dir     = '/home/ekaestne/PROJECTS/DATA/'; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf
dti_fsr_dir = '/home/mmilmcdRSI/data/proc_dti/';
fmr_fsr_dir = '/home/mmilmcdRSI/data/fsurf/';

prj_nme     = 'LanguageReorganization';
sbj_grp_nme = 'Lng_Reo_subject_list.csv';

%% Types
fib_typ = { 'FA' 'MD' };
scl_fct = [ 1    1000 ];

fib_cde = [101 102 103 104 105 106 107 108 109 110 115 116 117 118 119 120 121 122 123 133 134 135 136 137 138 141 142 143 144 145 146 147 148 149 150 1014 1024 2000 2001 2002 2003 2004];

%% Put together
sbj_grp = mmil_readtext([ prj_dir '/' 'SUBJECTS' '/' 'tmp' '/' prj_nme '/' sbj_grp_nme ]);
sbj_fsr = mmil_readtext([ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_freesurfer_names.csv' ]);
sbj_bld = mmil_readtext([ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_bold_names.csv' ]);

%%
for iFT = 1:numel(fib_typ)
    
    for iS = 1:size(sbj_grp,1) %
        
        fprintf( [ 'Subject #' num2str(iS) ' of ' num2str(size(sbj_grp,1)) ' : '  sbj_grp{iS,1} '\n' ] );
        
        sbj_fsr_ind = find( strcmpi( sbj_fsr(:,2) , sbj_grp{iS,1} ) );
        
        if ~isempty( sbj_fsr{sbj_fsr_ind,3} )
            
            fcfg = [];
            
            fcfg.ovr_dir = dti_fsr_dir;
            
            fcfg.sbj_nme     = sbj_grp{iS};
            fcfg.sbj_fsr_dir = sbj_fsr{iS,3};
            
            fcfg.atl_dir = 'AtlasTrack';
            fcfg.map_dir = 'fiber_maps';
            fcfg.fib_typ = fib_typ{iFT};
            fcfg.scl_fct = scl_fct(iFT);
            
            fcfg.thr_prb = 0.08;
            
            fcfg.fib_cde = fib_cde;
            
            fcfg.min_val = 1e-06;
            
            [ fib_dta_hld , tot_lbl] = ejk_extract_fibers(fcfg);
            
            if ~isempty(fib_dta_hld); fib_dta(iS,:) = fib_dta_hld; fib_nme_hld = tot_lbl; else fib_dta(iS,:) = nan(1,size(fib_dta,2)); end
            
        else
            
            fib_dta(iS,:) = nan(1,size(fib_dta,2));
            
        end
    end
    
    sve_sbj_nme = sbj_fsr(:,2);
    
    ejk_chk_dir( [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'Fibers' '/' ] )
    save(     [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'Fibers' '/' 'Fibers' '_' mmil_spec_char(fib_typ{iFT},{'.' '-' ' '}) '.mat'],'sve_sbj_nme','fib_nme_hld','fib_dta');
    cell2csv( [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'Fibers' '/' 'Fibers' '_' mmil_spec_char(fib_typ{iFT},{'.' '-' ' '}) '.csv'],[ ['SubjId' ; sve_sbj_nme(1:size(fib_dta,1))] [fib_nme_hld ; num2cell(fib_dta)] ])
    
    clear fib_dta
    
end

toc

