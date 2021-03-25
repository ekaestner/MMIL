

function ejk_extract_grey_thickness_wrapper

prj_dir = '/home/ekaestne/PROJECTS/';

prc_nme = { '' '.a2009s' '.split' };

%%
sbj_nme = mmil_readtext( [ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_subject_name.csv' ] );
fsr_nme = mmil_readtext( [ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_freesurfer_recons.csv' ] );

%%
for iPN = 1:numel(prc_nme)
   
    fcfg = [];

    fcfg.prj_dir = prj_dir;
    
    fcfg.prc_nme = prc_nme{iPN};
    
    for iS = 1:size(sbj_nme,1)
        
        fcfg.sbj_nme     = fsr_nme{iS,1};
        fcfg.ovr_dir     = fsr_nme{iS,2};
        fcfg.sbj_fsr_dir = fsr_nme{iS,3};
        
        [ gry_thk_dta_hld , tot_lbl ]= ejk_extract_grey_thickness(fcfg);
        if ~isempty(tot_lbl); tot_lbl = cellfun( @(x) mmil_spec_char( x , {'-'} ) , tot_lbl , 'uni' , 0 ); end
        
        try
            if ~isempty(gry_thk_dta_hld)
                gry_thk_dta(iS,:) = gry_thk_dta_hld;
                fprintf( [ sbj_nme{iS,1} ' : ' 'Data Loaded\n' ] )
            else
                gry_thk_dta(iS,:) = nan(1,size(gry_thk_dta,2));
                fprintf( [ sbj_nme{iS,1} ' : ' 'Missing\n' ] )
            end
        catch
            gry_thk_dta(iS,:) = nan(1,size(gry_thk_dta,2));
            fprintf( [ sbj_nme{iS,1} ' : ' 'Missing\n' ] )
        end
        
    end
    
    %% Save data
    sve_sbj_nme = sbj_nme(:,1);
    cell2csv([prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'MRI' '_' 'thickness' '_' 'aparc' '_' mmil_spec_char(fcfg.prc_nme,{'.'}) '.csv'],[ ['SubjId' ; sve_sbj_nme(1:size(gry_thk_dta,1))] [tot_lbl ; num2cell(gry_thk_dta)] ]);
    clear gry_thk_dta tot_lbl
    
end

end