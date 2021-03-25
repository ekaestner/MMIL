function ejk_project_pull_roi(cfg)

%% Load Data
dta_hld = mmil_readtext( [ cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' cfg.dta_nme '_' 'aparc' '_' mmil_spec_char(cfg.prc_nme,{'.'}) '.csv' ] );

%% Find Data
roi_hld = [cell( size(cfg.sbj_nme,1) , 1 ) num2cell(nan( size(cfg.sbj_nme,1) , size(dta_hld,2)-1 )) ];

for iR = 1:size( roi_hld , 1 )
    
    roi_hld{iR,1} = cfg.sbj_nme{iR,1};
    
    row_loc = find( strcmpi( dta_hld(:,1) , cfg.sbj_nme{iR,1} ) );
    roi_hld( iR , 2:end ) = dta_hld( row_loc , 2:end );
    
end

%% Output
cell2csv( [ cfg.out_dir  '/' cfg.dta_nme '_' 'aparc' '_' mmil_spec_char(cfg.prc_nme,{'.'}) '_' cfg.prj_nme '.csv'] , [ dta_hld( 1 , : ) ; roi_hld ] )

end