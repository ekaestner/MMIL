load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

dta_dir = [ prj_dir '/' prj_nme '/' 'Data' ];

out_dir = [ prj_dir '/' prj_nme '/' 'SpecificCor' '/' 'LogisticRegression' '/' ]; ejk_chk_dir( out_dir );

cat_cut = -1.5;

%% Load Data
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
% fcfg.all_num = 1;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'subcort_vol_ICV_cor_QC.csv'];
fcfg.dta_col = 5;
% fcfg.all_num = 1;
[ mri_dta, mri_dta_sbj, mri_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'fiber_FA_QC.csv'];
fcfg.dta_col = 5;
% fcfg.all_num = 1;
[ fib_dta, fib_dta_sbj, fib_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'wmparc_FA_wm_aparc_annot_QC.csv'];
fcfg.dta_col = 5;
% fcfg.all_num = 1;
[ wmp_dta, wmp_dta_sbj, wmp_dta_col] = ejk_dta_frm( fcfg );

prd_dta_sbj = [ cog_dta_sbj ];

prd_dta_col = [ cog_dta_col cln_dta_col mri_dta_col fib_dta_col wmp_dta_col ];

prd_dta     = [ cog_dta     cln_dta     mri_dta     fib_dta     wmp_dta ];

%% 
bnt_pre_scr_col = find(strcmpi( prd_dta_col, 'bnt_raw_scr' ));
ant_pre_scr_col = find(strcmpi( prd_dta_col, 'ant_mem_raw_scr' ));

lft_hip_col = find(strcmpi( prd_dta_col, 'Left_Hippocampus' ));
rgh_hip_col = find(strcmpi( prd_dta_col, 'Right_Hippocampus' ));

lft_ifo_col = find(strcmpi( prd_dta_col, 'L_IFO' ));
rgh_ifo_col = find(strcmpi( prd_dta_col, 'R_IFO' ));

lft_ilf_col = find(strcmpi( prd_dta_col, 'L_ILF' ));
rgh_ilf_col = find(strcmpi( prd_dta_col, 'R_ILF' ));

lft_fus_col = find(strcmpi( prd_dta_col, 'lh_fusiform' ));
rgh_fus_col = find(strcmpi( prd_dta_col, 'rh_fusiform' ));

cog_nme = { 'BNT'           'ANT' };
cog_col = [ bnt_pre_scr_col ant_pre_scr_col ];

mes_nme = { 'Hippocampus' 'ILF'       'IFOF'      'Fusiform' };
mes_col = [ lft_hip_col   lft_ilf_col lft_ifo_col  lft_fus_col ];

%%
grp_nme = { 'tle_controls_pre_3T_allSurg_all' 'controls_pre_3T_allSurg_all' 'tle_pre_3T_allSurg_left' 'tle_pre_3T_allSurg_right' };
grp_col = { rgb('black')                      rgb('light grey')              rgb('royal purple')       rgb('light magenta') };

grp_num = [ grp.controls_pre_3T_allSurg_all ; grp.tle_pre_3T_allSurg_left ; grp.tle_pre_3T_allSurg_right ];


%% 
for iT = 1:numel(cog_nme)
    for iN = 1:numel(mes_nme)
              
        % Data Gather
        for iG = 1:numel(grp_nme)
            if iG==1
                use_ind = intersect(grp.(grp_nme{iG}),grp_num);
            else
                use_ind = grp.(grp_nme{iG});
            end
            ydt{iG} = cell2mat(prd_dta( use_ind, mes_col(iN) ));
            xdt{iG} = cell2mat(prd_dta( use_ind, cog_col(iT) ));
        end
        
        % Data Plot
        fcfg = [];
        
        fcfg.xdt     = xdt;
        fcfg.ydt     = ydt;
                
        fcfg.fce_col = grp_col;
        fcfg.edg_col = { rgb('black') rgb('black') rgb('black') rgb('black') };
                
        fcfg.ylb = { mes_nme{iN}  };
        fcfg.xlb = { cog_nme{iT} };
        
        fcfg.trd_lne = ones(1,numel(grp_nme));
        
        fcfg.out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/atl_nme/figures/Figure3/';
        fcfg.out_nme = [ cog_nme{iT} '_' mes_nme{iN} ];
        
        ejk_scatter(fcfg)
        
        
    end    
end













