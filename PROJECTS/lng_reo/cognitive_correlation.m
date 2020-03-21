
lft_ind = find(strcmpi( grp_fle(:,3) , 'L' ));
rgh_ind = find(strcmpi( grp_fle(:,3) , 'R' ));
% lft_ind = 1:5; % find(strcmpi( , ))
con_ind = find(strcmpi( grp_fle(:,3) , 'HC' ));

%% DATA LOAD
% fMRI %%%%%%%%%%%%%%%
fmr_dta = mmil_readtext( [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'fMRI' '/' 'fMRI_aparc_xa2009s_N_FF_xnzvoxels.csv' ] );
fmr_lbl = fmr_dta(1,2:end);
fmr_dta = cell2mat(fmr_dta(2:end,2:end));

% WMParc MD %%%%%%%%%%%%%%%
wmp_dta = mmil_readtext( [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'WMParc' '/' 'WMParc_aparc__MD.csv' ] );
wmp_lbl = wmp_dta(1,2:end);
wmp_dta = cell2mat(wmp_dta(2:end,2:end));

% Grey Thickness %%%%%%%%%%%%%%%
gry_thk_dta = mmil_readtext( [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'GreyThick' '/' 'GreyThick_aparc_xa2009s.csv' ] );
gry_thk_lbl = gry_thk_dta(1,2:end);
gry_thk_dta = cell2mat(gry_thk_dta(2:end,2:end));

% Fibers %%%%%%%%%%%%%%%
fib_tfa_dta = mmil_readtext( [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'Fibers' '/' 'Fibers_FA.csv' ] );
fib_tfa_lbl = fib_tfa_dta(1,2:end);
fib_tfa_dta = cell2mat(fib_tfa_dta(2:end,2:end));

fib_tfa_lbl = fib_tfa_lbl(1:37);
fib_tfa_dta = fib_tfa_dta(:,1:37);

%% Laterality
% fMRI %%%%%%%%%%%%%%%
[ fmr_lat_dta , fmr_lat_lbl ]         = ejk_create_laterality_index( fmr_dta , fmr_lbl );

% WMParc MD %%%%%%%%%%%%%%%
[ wmp_lat_dta , wmp_lat_lbl ]         = ejk_create_laterality_index( wmp_dta , wmp_lbl );

% Grey Thickness %%%%%%%%%%%%%%%
[ gry_thk_lat_dta , gry_thk_lat_lbl ] = ejk_create_laterality_index( gry_thk_dta , gry_thk_lbl );

% Fibers %%%%%%%%%%%%%%%
[ fib_tfa_lat_dta , fib_tfa_lat_lbl ] = ejk_create_laterality_index( fib_tfa_dta , fib_tfa_lbl );

%% Z-Scores

%% Cognitive Put Together
sbj_cog_use.log_mem_nor_scr_one = sbj_cog.log_mem_nor_scr_one;
sbj_cog_use.log_mem_nor_scr_two = sbj_cog.log_mem_nor_scr_two;
sbj_cog_use.cvl_lfr_nor_scr_pst = sbj_cog.cvl_lfr_nor_scr_pst;
sbj_cog_use.vp1_nor_scr         = sbj_cog.vp1_nor_scr;
sbj_cog_use.vp2_nor_scr         = sbj_cog.vp2_nor_scr;
sbj_cog_use.bnt_nor_scr         = sbj_cog.bnt_nor_scr;
sbj_cog_use.ant_mem_raw_scr     = sbj_cog.ant_mem_raw_scr;
sbj_cog_use.cat_flu_nor_scr     = sbj_cog.cat_flu_nor_scr;
sbj_cog_use.cvl_tot_nor_scr     = sbj_cog.cvl_tot_nor_scr;
sbj_cog_use.ltr_tot_nor_scr     = sbj_cog.ltr_tot_nor_scr;
sbj_cog_use.swt_cor_nor_scr     = sbj_cog.swt_cor_nor_scr;
sbj_cog_use.swt_acc_nor_scr     = sbj_cog.swt_acc_nor_scr;

%% Put Together Neuroimaging Data
% fMRI %%%%%%%%%%%%%%%
fmr_cor_dta.sbj_nme = grp_fle(:,1);
for iC = 1:numel(fmr_lbl)     
    fmr_cor_dta.(mmil_spec_char(fmr_lbl{iC},{'-'})) = fmr_dta(:,iC);
end

fmr_lat_cor_dta.sbj_nme = grp_fle(:,1);
for iC = 1:numel(fmr_lat_lbl)     
    fmr_lat_cor_dta.(mmil_spec_char(fmr_lat_lbl{iC},{'-'})) = fmr_lat_dta(:,iC);
end

% WMParc MD %%%%%%%%%%%%%%%
wmp_cor_dta.sbj_nme = grp_fle(:,1);
for iC = 1:numel(wmp_lbl)     
    wmp_cor_dta.(mmil_spec_char(wmp_lbl{iC},{'-'})) = wmp_dta(:,iC);
end

wmp_lat_cor_dta.sbj_nme = grp_fle(:,1);
for iC = 1:numel(wmp_lat_lbl)     
    wmp_lat_cor_dta.(mmil_spec_char(wmp_lat_lbl{iC},{'-'})) = wmp_lat_dta(:,iC);
end

% Grey Thickness %%%%%%%%%%%%%%%
gry_thk_cor_dta.sbj_nme = grp_fle(:,1);
for iC = 1:numel(gry_thk_lbl)     
    gry_thk_cor_dta.(mmil_spec_char(gry_thk_lbl{iC},{'-'})) = gry_thk_dta(:,iC);
end

gry_thk_lat_cor_dta.sbj_nme = grp_fle(:,1);
for iC = 1:numel(gry_thk_lat_lbl)     
    gry_thk_lat_cor_dta.(mmil_spec_char(gry_thk_lat_lbl{iC},{'-'})) = gry_thk_lat_dta(:,iC);
end

% Fibers %%%%%%%%%%%%%%%
fib_tfa_cor_dta.sbj_nme = grp_fle(:,1);
for iC = 1:numel(fib_tfa_lbl)     
    fib_tfa_cor_dta.(mmil_spec_char(fib_tfa_lbl{iC},{'-'})) = fib_tfa_dta(:,iC);
end

fib_tfa_lat_cor_dta.sbj_nme = grp_fle(:,1);
for iC = 1:numel(fib_tfa_lat_lbl)     
    fib_tfa_lat_cor_dta.(mmil_spec_char(fib_tfa_lat_lbl{iC},{'-'})) = fib_tfa_lat_dta(:,iC);
end

%% Neuroimaging by Neuroimaging
fcfg = [];

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Correlations' '/' 'fMRICorrelation' '/'];

fcfg.dta_lbl = { 'fmr_by_tfa_dta' 'fmr_by_wmd_dta' 'fmr_by_gry_dta' };
fcfg.xdt     = { fib_tfa_cor_dta  wmp_cor_dta      gry_thk_cor_dta  }; % Secondary Folders
fcfg.ydt     = { fmr_cor_dta      fmr_cor_dta      fmr_cor_dta      }; % Master Folders

fcfg.grp     = grp_fle;
fcfg.grp_col = [ 1                    2 ];
fcfg.grp_nme = { { 'EPD'    'HC' }    { 'L'    'R'   } };
fcfg.grp_clr = { { 'orange' 'black' } { 'blue' 'red' } };

ejk_roi_scatter(fcfg)

%% Cognitive by Neuroimaging
fcfg = [];

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Correlations' '/' 'CognitiveCorrelation' ];

fcfg.dta_lbl = { 'wmp_cor_dta' 'gry_thk_cor_dta' 'fib_tfa_cor_dta'       }; % 'fmr_cor_dta'
fcfg.xdt     = { sbj_cog_use   sbj_cog_use       sbj_cog_use     }; % Secondary Folders % sbj_cog_use
fcfg.ydt     = { wmp_cor_dta   gry_thk_cor_dta   fib_tfa_cor_dta }; % Master Folders % fmr_cor_dta

fcfg.grp     = grp_fle;
fcfg.grp_col = [ 1                    2 ];
fcfg.grp_nme = { { 'EPD'    'HC' }    { 'L'    'R'   } };
fcfg.grp_clr = { { 'orange' 'black' } { 'blue' 'red' } };

ejk_roi_scatter(fcfg)

%% Neuroimaging Laterality by Neuroimaging Laterality
fcfg = [];

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Correlations' '/' 'fMRICorrelation' '/'];

fcfg.dta_lbl = { 'fmr_by_tfa_dta'    'fmr_by_wmd_dta' 'fmr_by_gry_dta' };
fcfg.xdt     = { fib_tfa_lat_cor_dta wmp_lat_cor_dta  gry_thk_lat_cor_dta  }; % Secondary Folders
fcfg.ydt     = { fmr_lat_cor_dta     fmr_lat_cor_dta  fmr_lat_cor_dta      }; % Master Folders

fcfg.grp     = grp_fle;
fcfg.grp_col = [ 1                    2 ];
fcfg.grp_nme = { { 'EPD'    'HC' }    { 'L'    'R'   } };
fcfg.grp_clr = { { 'orange' 'black' } { 'blue' 'red' } };

ejk_roi_scatter(fcfg)

%% Cognitive by Neuroimaging Laterality
fcfg = [];

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Correlations' '/' 'fMRICorrelation_Laterality' ];

fcfg.dta_lbl = { 'fmr_cor_dta'   'wmp_cor_dta'     'gry_thk_cor_dta'   'fib_tfa_cor_dta'   }; % 
fcfg.xdt     = { sbj_cog_use     sbj_cog_use       sbj_cog_use         sbj_cog_use         }; % Secondary Folders % 
fcfg.ydt     = { fmr_lat_cor_dta wmp_lat_cor_dta   gry_thk_lat_cor_dta fib_tfa_lat_cor_dta }; % Master Folders

fcfg.grp     = grp_fle;
fcfg.grp_col = [ 1                    2 ];
fcfg.grp_nme = { { 'EPD'    'HC' }    { 'L'    'R'   } };
fcfg.grp_clr = { { 'orange' 'black' } { 'blue' 'red' } };

ejk_roi_scatter(fcfg)