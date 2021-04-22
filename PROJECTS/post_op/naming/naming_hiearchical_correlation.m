
load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

run_grp = { 'tle_post_3T_ATLonly_left' };

dta_dir = [ prj_dir '/' prj_nme '/' 'Data' ];

out_dir = [ prj_dir '/' prj_nme '/' 'SpecificCor' '/' 'HiearchicalRegression' '/' ]; ejk_chk_dir( out_dir );

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

prd_dta_col = [ cog_dta_col cln_dta_col mri_dta_col fib_dta_col wmp_dta_col ];

prd_dta     = [ cog_dta     cln_dta     mri_dta     fib_dta     wmp_dta ];

%% Setup
out_nme     = { 'BNT'              'ANT' };

dep_var_nme = { 'xbnt_raw_scr_pst' 'xant_mem_raw_scr_pst' };
dep_var_dta = prd_dta( grp.tle_post_3T_ATLonly_left, find(ismember(prd_dta_col,dep_var_nme)));

prd_var_nme{1}{1} = sort({ 'xbnt_raw_scr' 'xLeft_Hippocampus' });
prd_var_nme{1}{2} = sort({ 'xrh_fusiform' 'xL_ILF' });
prd_var_dta{1}    = [ prd_dta( grp.tle_post_3T_ATLonly_left, find(ismember(prd_dta_col,prd_var_nme{1}{1}))) ...    
                      prd_dta( grp.tle_post_3T_ATLonly_left, find(ismember(prd_dta_col,prd_var_nme{1}{2}))) ];

prd_var_nme{2}{1} = sort({ 'xant_mem_raw_scr' 'xLeft_Hippocampus' });
prd_var_nme{2}{2} = sort({ 'xrh_fusiform' 'xL_ILF' });
prd_var_dta{2} = [ prd_dta( grp.tle_post_3T_ATLonly_left, find(ismember(prd_dta_col,prd_var_nme{2}{1}))) ...
                   prd_dta( grp.tle_post_3T_ATLonly_left, find(ismember(prd_dta_col,prd_var_nme{2}{2}))) ];

for iO = 1:numel(out_nme)

    dep_var.sbj_nme = cog_dta_sbj(grp.tle_post_3T_ATLonly_left);
    dep_var.(dep_var_nme{iO}) = dep_var_dta(:,iO);
    
    prd_var.sbj_nme = cog_dta_sbj(grp.tle_post_3T_ATLonly_left);
    nme_hld = cat(2,prd_var_nme{iO}{:});
    for iD = 1:numel(nme_hld)
        prd_var.(nme_hld{iD}) = prd_var_dta{iO}(:,iD);
    end
    for iG = 1:numel(prd_var_nme{iO})        
        prd_nme.(['blk' '_' num2str(iG)]) = prd_var_nme{iO}{iG};
    end
   
    ejk_chk_dir(out_dir)
    save( [out_dir '/' 'dep_var.mat' ], 'dep_var' )
    save( [out_dir '/' 'prd_var.mat' ], 'prd_var' )
    save( [out_dir '/' 'prd_nme.mat' ], 'prd_nme' )
    
    clear dep_var prd_var prd_nme
    
end

%% Save Data
cell2csv( '/home/ekaestne/Downloads/postoperative_language_subjects_ant.csv',[ [ cog_dta_sbj(grp.tle_post_3T_ATLonly_left,1) ; cog_dta_sbj(grp.tle_post_3T_ATLonly_left,1) ] ...
  [ cog_dta(grp.tle_post_3T_ATLonly_left,2)     ; cog_dta(grp.tle_post_3T_ATLonly_left,2) ]  ...
  [ cog_dta(grp.tle_post_3T_ATLonly_left,5)     ; cog_dta(grp.tle_post_3T_ATLonly_left,5) ]  ]);








