clear; clc;

%
prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'scratch';

%
frs_dta = '/home/mmilmcd/data/FSRECONS/';

%
red_fle = 'sbj000_total_2019_03_27.csv';

%% Redcap Data
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.red_fle = red_fle;
[sbj_dem , sbj_sze , sbj_scn , sbj_cog] = mmil_load_redcap(fcfg);

%% Freesurfer names
frs_nme = mmil_readtext( [ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_freesurfer_names.csv' ] );

%%
for iS = 169:numel( sbj_dem.sbj_nme )

    ejk_chk_dir( [ prj_dir '/' 'DATA' '/' sbj_dem.sbj_nme{iS} ] )
    ejk_chk_dir( [ prj_dir '/' 'DATA' '/' sbj_dem.sbj_nme{iS} '/' 'parcellations'] )
    
    if ~isempty(frs_nme{iS,4})
        fprintf(['Copying ' sbj_dem.sbj_nme{iS} ' : ' frs_nme{iS,4} '\n'])
        copyfile( [ frs_dta '/' frs_nme{iS,4} '/' 'label' '/' 'lh.aparc.a2009s.annot' ] , [ prj_dir '/' 'DATA' '/' sbj_dem.sbj_nme{iS} '/' 'parcellations' '/' 'lh.aparc.a2009s.annot' ] );
        copyfile( [ frs_dta '/' frs_nme{iS,4} '/' 'label' '/' 'rh.aparc.a2009s.annot' ] , [ prj_dir '/' 'DATA' '/' sbj_dem.sbj_nme{iS} '/' 'parcellations' '/' 'rh.aparc.a2009s.annot' ] );
        copyfile( [ frs_dta '/' frs_nme{iS,4} '/' 'label' '/' 'lh.aparc.annot' ]        , [ prj_dir '/' 'DATA' '/' sbj_dem.sbj_nme{iS} '/' 'parcellations' '/' 'lh.aparc.annot' ] );
        copyfile( [ frs_dta '/' frs_nme{iS,4} '/' 'label' '/' 'rh.aparc.annot' ]        , [ prj_dir '/' 'DATA' '/' sbj_dem.sbj_nme{iS} '/' 'parcellations' '/' 'rh.aparc.annot' ] );
    end

end
    