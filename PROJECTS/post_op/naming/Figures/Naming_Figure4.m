load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

dta_dir = [ prj_dir '/' prj_nme '/' 'Data' ];

out_dir = [ prj_dir '/' prj_nme '/' 'Manuscript' '/' 'Figures' '/' 'Figure2' '/' ]; ejk_chk_dir( out_dir );

cog_dta_nme = [ dta_dir '/' 'Cognitive'                '_' 'QC' '.csv'];

%% Load Data
load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

dta_dir = [ prj_dir '/' prj_nme '/' 'Data' ];

out_dir = [ prj_dir '/' prj_nme '/' 'Manuscript' '/' 'Figures' '/' 'Figure3' '/' ]; ejk_chk_dir( out_dir );

cog_dta_nme = [ dta_dir '/' 'Cognitive'                '_' 'QC' '.csv'];
cln_dta_nme = [ dta_dir '/' 'Clinical'                          '.csv'];
mri_dta_nme = [ dta_dir '/' 'subcort_vol_ICV_cor'      '_' 'QC' '.csv'];
fib_dta_nme = [ dta_dir '/' 'fiber_FA'                 '_' 'QC' '.csv'];
wmp_dta_nme = [ dta_dir '/' 'wmparc_FA_wm_aparc_annot' '_' 'QC' '.csv'];

%% Load Data
cog_dta = mmil_readtext(cog_dta_nme);
cog_dta_col = ejk_fix_column_names(cog_dta(1,2:end));
cog_dta_sbj = cog_dta(2:end,1);
cog_dta     = cog_dta(2:end,2:end);

cln_dta = mmil_readtext(cln_dta_nme);
cln_dta_col = ejk_fix_column_names(cln_dta(1,2:end));
cln_dta_sbj = cln_dta(2:end,1);
cln_dta     = cln_dta(2:end,2:end);

mri_dta = mmil_readtext(mri_dta_nme);
mri_dta_col = ejk_fix_column_names(mri_dta(1,5:end));
mri_dta_sbj = mri_dta(2:end,1);
mri_dta     = mri_dta(2:end,5:end);

fib_dta = mmil_readtext(fib_dta_nme);
fib_dta_col = ejk_fix_column_names(fib_dta(1,5:end));
fib_dta_sbj = fib_dta(2:end,1);
fib_dta     = fib_dta(2:end,5:end);

wmp_dta = mmil_readtext(wmp_dta_nme);
wmp_dta_col = ejk_fix_column_names(wmp_dta(1,5:end));
wmp_dta_sbj = wmp_dta(2:end,1);
wmp_dta     = wmp_dta(2:end,5:end);

prd_dta_sbj = [ cog_dta_sbj ];

prd_dta_col = [ cog_dta_col cln_dta_col mri_dta_col fib_dta_col wmp_dta_col ];

prd_dta     = [ cog_dta     cln_dta     mri_dta     fib_dta     wmp_dta ];

%% 
bnt_pst_scr_col = find(strcmpi( prd_dta_col, 'xbnt_raw_scr_pst' ));
ant_pst_scr_col = find(strcmpi( prd_dta_col, 'xant_mem_raw_scr_pst' ));

lft_hip_col = find(strcmpi( prd_dta_col, 'xLeft_Hippocampus' ));
rgh_hip_col = find(strcmpi( prd_dta_col, 'xRight_Hippocampus' ));

lft_ifo_col = find(strcmpi( prd_dta_col, 'xL_IFO' ));
rgh_ifo_col = find(strcmpi( prd_dta_col, 'xR_IFO' ));

lft_ilf_col = find(strcmpi( prd_dta_col, 'xL_ILF' ));
rgh_ilf_col = find(strcmpi( prd_dta_col, 'xR_ILF' ));

lft_fus_col = find(strcmpi( prd_dta_col, 'xlh_fusiform' ));
rgh_fus_col = find(strcmpi( prd_dta_col, 'xrh_fusiform' ));

cog_nme = { 'BNT'           'ANT' };
cog_col = [ bnt_pst_scr_col ant_pst_scr_col ];

mes_nme = { 'Hippocampus' 'ILF'       'IFOF'      'Fusiform' };
mes_col = [ lft_hip_col   lft_ilf_col lft_ifo_col rgh_fus_col ];

%%
grp_nme = { 'tle_post_3T_ATLonly_left' 'tle_post_3T_ATLonly_right' };
grp_col = { rgb('royal purple')       rgb('light magenta') };

%% 
for iT = 1:numel(cog_nme)
    for iN = 1:numel(mes_nme)
       
        % Data Gather
        for iG = 1:numel(grp_nme)
           ydt{iG} = cell2mat(prd_dta( grp.(grp_nme{iG}), mes_col(iN) )); 
           xdt{iG} = cell2mat(prd_dta( grp.(grp_nme{iG}), cog_col(iT) )); 
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
        
        fcfg.out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Naming/Manuscript/Figures/Figure4/';
        fcfg.out_nme = [ cog_nme{iT} '_' mes_nme{iN} ];
        
        ejk_scatter(fcfg)
        
        
    end    
end






