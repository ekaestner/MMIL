
% %%%%%%%%%%%%%%%
wmp_dta = mmil_readtext( [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'WMParc' '/' 'WMParc_aparc__MD.csv' ] );
wmp_lbl = wmp_dta(1,2:end);
wmp_dta = cell2mat(wmp_dta(2:end,2:end));

wmp_dta_tot = nanmean( wmp_dta , 2 );

% %%%%%%%%%%%%%%%
lft_ind = find(strcmpi( grp_fle(:,3) , 'L' ));
rgh_ind = find(strcmpi( grp_fle(:,3) , 'R' ));
% lft_ind = 1:5; % find(strcmpi( , ))
con_ind = find(strcmpi( grp_fle(:,3) , 'HC' ));

% Z-Score %%%%%%%%%%%%%%%
fcfg  = [];
fcfg.con_ind = con_ind;
wmp_dta_zsc = ejk_create_zscore(fcfg,wmp_dta);

fcfg  = [];
fcfg.con_ind = con_ind;
wmp_dta_tot_zsc = ejk_create_zscore(fcfg,wmp_dta_tot);

%% Total Activation
% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw Data 
lft_tot = wmp_dta_tot( lft_ind , : );
rgh_tot = wmp_dta_tot( rgh_ind , : );
con_tot = wmp_dta_tot( con_ind , : ) ;

fcfg = [];

fcfg.xdt     = { 1            2            3 };
fcfg.ydt     = { lft_tot      rgh_tot      con_tot };

% fcfg.ylm     = [ 0 3000 ];

fcfg.fce_col = { rgb('blue')  rgb('red')   rgb('greenish grey') };
fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };

fcfg.xlb = '';
fcfg.ylb = 'Total wmpI Activation';

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'WMParcMD' '/' 'Desikan' '/'];
fcfg.out_nme = 'Desikan_OverallActivation';

ejk_scatter(fcfg)

% Z-Score
lft_tot_zsc = wmp_dta_tot_zsc( lft_ind , : );
rgh_tot_zsc = wmp_dta_tot_zsc( rgh_ind , : );
con_tot_zsc = wmp_dta_tot_zsc( con_ind , : );

fcfg = [];

fcfg.xdt     = { 1            2            3 };
fcfg.ydt     = { lft_tot_zsc  rgh_tot_zsc  con_tot_zsc };

fcfg.fce_col = { rgb('blue')  rgb('red')   rgb('greenish grey') };
fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };

fcfg.xlb = '';
fcfg.ylb = 'Total wmpI Zscore';

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'WMParcMD' '/' 'Desikan' '/'];
fcfg.out_nme = 'Desikan_OverallActivation_ZScore';

ejk_scatter(fcfg)

% Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw Data 
fcfg = [];

fcfg.ind     = { lft_ind rgh_ind con_ind };
fcfg.ind_nme = { 'left'  'right' 'HC' };

fcfg.dta_lbl = {'Total Activation'};

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'WMParcMD' '/' 'Desikan' '/'];
fcfg.out_nme = 'Desikan_OverallActivation';

ejk_grp_stt( fcfg , wmp_dta_tot )

%% ROIs
% Raw Data
for iRO = 1:size(wmp_lbl,2)
    
    %
    lft_roi = wmp_dta( lft_ind , iRO );
    rgh_roi = wmp_dta( rgh_ind , iRO );
    con_roi = wmp_dta( con_ind , iRO );
    
    %
    fcfg = [];
    
    fcfg.xdt     = { 1            2            3 };
    fcfg.ydt     = { lft_roi      rgh_roi      con_roi };
    
    fcfg.fce_col = { rgb('blue')  rgb('red')   rgb('greenish grey') };
    fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };
    
    fcfg.xlb = '';
    fcfg.ylb = 'wmpI Activation';
    
    fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'WMParcMD' '/' 'Desikan' '/' 'Activation'];
    fcfg.out_nme = mmil_spec_char( wmp_lbl{iRO} , {'-'} );
    
    ejk_scatter(fcfg)
    
end

% Z-Score
for iRO = 1:size(wmp_lbl,2)
    
    %
    lft_roi = wmp_dta_zsc( lft_ind , iRO );
    rgh_roi = wmp_dta_zsc( rgh_ind , iRO );
    con_roi = wmp_dta_zsc( con_ind , iRO );
    
    %
    if ~ (all(isnan(lft_roi)) && all(isnan(rgh_roi)) && all(isnan(con_roi)))
        
        fcfg = [];
        
        fcfg.xdt     = { 1            2            3 };
        fcfg.ydt     = { lft_roi      rgh_roi      con_roi };
        
        fcfg.fce_col = { rgb('blue')  rgb('red')   rgb('greenish grey') };
        fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };
        
        fcfg.xlb = '';
        fcfg.ylb = 'wmpI Zscore';
        
        fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'WMParcMD' '/' 'Desikan' '/' 'Activation_zscore'];
        fcfg.out_nme = [ mmil_spec_char( wmp_lbl{iRO} , {'-'} ) '_' 'zscore'];
        
        ejk_scatter(fcfg)
        
    end
    
end

% Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw Data 
fcfg = [];

fcfg.ind     = { lft_ind rgh_ind con_ind };
fcfg.ind_nme = { 'left'  'right' 'HC' };

fcfg.dta_lbl = wmp_lbl(2:end)';

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'WMParcMD' '/' 'Desikan' '/'];
fcfg.out_nme = 'Desikan_ROIActivation';

ejk_grp_stt( fcfg , wmp_dta(:,2:end) )




