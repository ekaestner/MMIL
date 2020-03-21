%%
% %%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%
wmp_dta = mmil_readtext( [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'WMParc' '/' 'WMParc_aparc__MD.csv' ] );
wmp_lbl = wmp_dta(1,2:end);
wmp_dta = cell2mat(wmp_dta(2:end,2:end));

% %%%%%%%%%%%%%%%
[ wmp_lat_dta , wmp_lat_lbl ] = ejk_create_laterality_index( wmp_dta , wmp_lbl );

tot_dta(:,1) = nanmean( wmp_dta( : , string_find( wmp_lbl , {'lh_'} ) ) , 2 );
tot_dta(:,2) = nanmean( wmp_dta( : , string_find( wmp_lbl , {'rh_'} ) ) , 2 );
tot_lbl = { 'lh_dta' 'rh_dta' };
[ wmp_tot_lat_dta , wmp_tot_lat_lbl ] = ejk_create_laterality_index( tot_dta , tot_lbl );

% %%%%%%%%%%%%%%%
lft_ind = find(strcmpi( grp_fle(:,3) , 'L' ));
rgh_ind = find(strcmpi( grp_fle(:,3) , 'R' ));
% lft_ind = 1:5; % find(strcmpi( , ))
con_ind = find(strcmpi( grp_fle(:,3) , 'HC' ));

% Z-Score %%%%%%%%%%%%%%%
fcfg  = [];
fcfg.con_ind = con_ind;
wmp_lat_dta_zsc = ejk_create_zscore(fcfg,wmp_lat_dta);

fcfg  = [];
fcfg.con_ind = con_ind;
wmp_tot_lat_dta_zsc = ejk_create_zscore(fcfg,wmp_tot_lat_dta);

%% Total Laterality
% Raw Data
lft_tot = wmp_tot_lat_dta( lft_ind , : );
rgh_tot = wmp_tot_lat_dta( rgh_ind , : );
con_tot = wmp_tot_lat_dta( con_ind , : );

fcfg = [];

fcfg.xdt     = { 1            2            3 };
fcfg.ydt     = { lft_tot      rgh_tot      con_tot };

fcfg.ylm     = [ -1.1 1.1 ];

fcfg.fce_col = { rgb('blue')  rgb('red')   rgb('greenish grey') };
fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };

fcfg.xlb = '';
fcfg.ylb = 'Total WMParcMD Laterality';

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'WMParcMD' '/' 'Desikan' '/' ];
fcfg.out_nme = 'Desikan_OverallLaterality';

ejk_scatter(fcfg)

% Z-Score
lft_tot_zsc = wmp_tot_lat_dta_zsc( lft_ind , : );
rgh_tot_zsc = wmp_tot_lat_dta_zsc( rgh_ind , : );
con_tot_zsc = wmp_tot_lat_dta_zsc( con_ind , : );

fcfg = [];

fcfg.xdt     = { 1            2            3 };
fcfg.ydt     = { lft_tot_zsc  rgh_tot_zsc  con_tot_zsc };

fcfg.fce_col = { rgb('blue')  rgb('red')   rgb('greenish grey') };
fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };

fcfg.xlb = '';
fcfg.ylb = 'Total WMParcMD Laterality Zscore';

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'WMParcMD' '/' 'Desikan' '/' ];
fcfg.out_nme = 'Desikan_OverallLaterality_ZScore';

ejk_scatter(fcfg)

% Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw Data 
fcfg = [];

fcfg.ind     = { lft_ind rgh_ind con_ind };
fcfg.ind_nme = { 'left'  'right' 'HC' };

fcfg.dta_lbl = {'Total Laterality'};

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'WMParcMD' '/' 'Desikan' '/'];
fcfg.out_nme = 'Desikan_OverallLaterality';

ejk_grp_stt( fcfg , wmp_tot_lat_dta )

%% ROIs
% Raw Data
for iRO = 1:size(wmp_lat_lbl,2)
    
    %
    lft_roi = wmp_lat_dta( lft_ind , iRO );
    rgh_roi = wmp_lat_dta( rgh_ind , iRO );
    con_roi = wmp_lat_dta( con_ind , iRO );
    
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
    fcfg.ylb = 'WMParcMD Laterality';
    
    fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'WMParcMD' '/' 'Desikan' '/' 'Laterality' '/' ];
    fcfg.out_nme = mmil_spec_char( wmp_lat_lbl{iRO} , {'-'} );
    
    ejk_scatter(fcfg)
    
end

% Z-Score
for iRO = 1:size(wmp_lat_lbl,2)
    
    %
    lft_roi = wmp_lat_dta_zsc( lft_ind , iRO );
    rgh_roi = wmp_lat_dta_zsc( rgh_ind , iRO );
    con_roi = wmp_lat_dta_zsc( con_ind , iRO );
    
    %
    if ~ (all(isnan(lft_roi)) && all(isnan(rgh_roi)) && all(isnan(con_roi)))
        
        fcfg = [];
        
        fcfg.xdt     = { 1            2            3 };
        fcfg.ydt     = { lft_roi      rgh_roi      con_roi };
        
        fcfg.fce_col = { rgb('blue')  rgb('red')   rgb('greenish grey') };
        fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };
        
        fcfg.xlb = '';
        fcfg.ylb = 'WMParcMD Laterality Zscore';
        
        fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'WMParcMD' '/' 'Desikan' '/' 'Laterality_zscore' '/' ];
        fcfg.out_nme = [ mmil_spec_char( wmp_lat_lbl{iRO} , {'-'} ) '_' 'zscore'];
        
        ejk_scatter(fcfg)
        
    end
    
end

% Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw Data 
fcfg = [];

fcfg.ind     = { lft_ind rgh_ind con_ind };
fcfg.ind_nme = { 'left'  'right' 'HC' };

fcfg.dta_lbl = wmp_lat_lbl';

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'WMParcMD' '/' 'Desikan' '/'];
fcfg.out_nme = 'Desikan_ROILaterality';

ejk_grp_stt( fcfg , wmp_lat_dta )



