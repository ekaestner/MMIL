%%
% %%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%
fib_tMD_dta = mmil_readtext( [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'Fibers' '/' 'Fibers_MD.csv' ] );
fib_tMD_lbl = fib_tMD_dta(1,2:end);
fib_tMD_dta = cell2mat(fib_tMD_dta(2:end,2:end));

% %%%%%%%%%%%%%%%
[ fib_tMD_lat_dta , fib_tMD_lat_lbl ] = ejk_create_laterality_index( fib_tMD_dta , fib_tMD_lbl );

tot_dta = fib_tMD_dta( : , 39:42 );
tot_lbl = fib_tMD_lbl( 1 , 39:42 );
[ fib_tMD_tot_lat_dta , fib_tMD_tot_lat_lbl ] = ejk_create_laterality_index( tot_dta , tot_lbl );

% %%%%%%%%%%%%%%%
lft_ind = find(strcmpi( grp_fle(:,3) , 'L' ));
rgh_ind = find(strcmpi( grp_fle(:,3) , 'R' ));
% lft_ind = 1:5; % find(strcmpi( , ))
con_ind = find(strcmpi( grp_fle(:,3) , 'HC' ));

% Z-Score %%%%%%%%%%%%%%%
fcfg  = [];
fcfg.con_ind = con_ind;
fib_tMD_lat_dta_zsc = ejk_create_zscore(fcfg,fib_tMD_lat_dta);

fcfg  = [];
fcfg.con_ind = con_ind;
fib_tMD_tot_lat_dta_zsc = ejk_create_zscore(fcfg,fib_tMD_tot_lat_dta);

%% Total Laterality
% Raw Data
for iRO = 1:size(fib_tMD_tot_lat_lbl,2)
    
    %
    lft_roi = fib_tMD_tot_lat_dta( lft_ind , iRO );
    rgh_roi = fib_tMD_tot_lat_dta( rgh_ind , iRO );
    con_roi = fib_tMD_tot_lat_dta( con_ind , iRO );
    
    lft_nan = sum(isnan( lft_roi ));
    rgh_nan = sum(isnan( rgh_roi ));
    con_nan = sum(isnan( con_roi ));
    
    %
    fcfg = [];
    
    fcfg.xdt     = { 1            2            3 };
    fcfg.ydt     = { lft_roi      rgh_roi      con_roi };
    
    fcfg.fce_col = { rgb('blue')  rgb('red')   rgb('greenish grey') };
    fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };
    
    fcfg.xlb = ['NaN: LTLE(' num2str(lft_nan) ') RTLE(' num2str(rgh_nan) ') HC(' num2str(con_nan) ')'];
    fcfg.ylb = 'fMRI Laterality';
    
    fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'fMRI' '/' 'Fibers_MD' '/' 'Laterality_total' '/' ];
    fcfg.out_nme = mmil_spec_char( fib_tMD_tot_lat_lbl{iRO} , {'-'} );
    
    ejk_scatter(fcfg)
    
end

% Z-Score
for iRO = 1:size(fib_tMD_tot_lat_lbl,2)
    
    %
    lft_roi = fib_tMD_tot_lat_dta_zsc( lft_ind , iRO );
    rgh_roi = fib_tMD_tot_lat_dta_zsc( rgh_ind , iRO );
    con_roi = fib_tMD_tot_lat_dta_zsc( con_ind , iRO );
    
    %
    if ~ (all(isnan(lft_roi)) && all(isnan(rgh_roi)) && all(isnan(con_roi)))
        
        fcfg = [];
        
        fcfg.xdt     = { 1            2            3 };
        fcfg.ydt     = { lft_roi      rgh_roi      con_roi };
        
        fcfg.fce_col = { rgb('blue')  rgb('red')   rgb('greenish grey') };
        fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };
        
        fcfg.xlb = '';
        fcfg.ylb = 'fMRI Laterality Zscore';
        
        fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'fMRI' '/' 'Fibers_MD' '/' 'Laterality_total_zscore' '/' ];
        fcfg.out_nme = [ mmil_spec_char( fib_tMD_tot_lat_lbl{iRO} , {'-'} ) '_' 'zscore'];
        
        ejk_scatter(fcfg)
        
    end
    
end

% Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw Data 
fcfg = [];

fcfg.ind     = { lft_ind rgh_ind con_ind };
fcfg.ind_nme = { 'left'  'right' 'HC' };

fcfg.dta_lbl = fib_tMD_tot_lat_lbl';

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'fMRI' '/' 'Fibers_MD' '/'];
fcfg.out_nme = 'Fibers_MD_OverallLaterality';

ejk_grp_stt( fcfg , fib_tMD_tot_lat_dta )

%% ROIs
% Raw Data
for iRO = 1:size(fib_tMD_lat_lbl,2)
    
    %
    lft_roi = fib_tMD_lat_dta( lft_ind , iRO );
    rgh_roi = fib_tMD_lat_dta( rgh_ind , iRO );
    con_roi = fib_tMD_lat_dta( con_ind , iRO );
    
    lft_nan = sum(isnan( lft_roi ));
    rgh_nan = sum(isnan( rgh_roi ));
    con_nan = sum(isnan( con_roi ));
    
    %
    fcfg = [];
    
    fcfg.xdt     = { 1            2            3 };
    fcfg.ydt     = { lft_roi      rgh_roi      con_roi };
    
    fcfg.fce_col = { rgb('blue')  rgb('red')   rgb('greenish grey') };
    fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };
    
    fcfg.xlb = ['NaN: LTLE(' num2str(lft_nan) ') RTLE(' num2str(rgh_nan) ') HC(' num2str(con_nan) ')'];
    fcfg.ylb = 'fMRI Laterality';
    
    fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'fMRI' '/' 'Fibers_MD' '/' 'Laterality' '/' ];
    fcfg.out_nme = mmil_spec_char( fib_tMD_lat_lbl{iRO} , {'-'} );
    
    ejk_scatter(fcfg)
    
end

% Z-Score
for iRO = 1:size(fib_tMD_lat_lbl,2)
    
    %
    lft_roi = fib_tMD_lat_dta_zsc( lft_ind , iRO );
    rgh_roi = fib_tMD_lat_dta_zsc( rgh_ind , iRO );
    con_roi = fib_tMD_lat_dta_zsc( con_ind , iRO );
    
    %
    if ~ (all(isnan(lft_roi)) && all(isnan(rgh_roi)) && all(isnan(con_roi)))
        
        fcfg = [];
        
        fcfg.xdt     = { 1            2            3 };
        fcfg.ydt     = { lft_roi      rgh_roi      con_roi };
        
        fcfg.fce_col = { rgb('blue')  rgb('red')   rgb('greenish grey') };
        fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };
        
        fcfg.xlb = '';
        fcfg.ylb = 'fMRI Laterality Zscore';
        
        fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'fMRI' '/' 'Fibers_MD' '/' 'Laterality_zscore' '/' ];
        fcfg.out_nme = [ mmil_spec_char( fib_tMD_lat_lbl{iRO} , {'-'} ) '_' 'zscore'];
        
        ejk_scatter(fcfg)
        
    end
    
end

% Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw Data 
fcfg = [];

fcfg.ind     = { lft_ind rgh_ind con_ind };
fcfg.ind_nme = { 'left'  'right' 'HC' };

fcfg.dta_lbl = fib_tMD_lat_lbl';

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'fMRI' '/' 'Fibers_MD' '/'];
fcfg.out_nme = 'Fibers_MD_ROILaterality';

ejk_grp_stt( fcfg , fib_tMD_lat_dta )



