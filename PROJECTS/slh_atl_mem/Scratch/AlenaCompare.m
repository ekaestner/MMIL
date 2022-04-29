out_put = [ prj_dir '/' prj_nme '/' 'InitialAnalysis'];

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/slh_atl_mem/Misc/AlenaCompare' '/' 'FinalDataset_ATLMemory.csv'];
fcfg.dta_col = 2;
[ aln_dta, aln_dta_sbj, aln_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/Data' '/' 'fiber_FA_238.csv'];
fcfg.dta_col = 2;
[ fib_dta, fib_dta_sbj, fib_dta_col] = ejk_dta_frm( fcfg );

%% Compare Subjects
cln_dta_sbj( grp.tle_controls_pre_allT_allSurg_all )

ejk_hld = cln_dta_sbj( grp.tle_post_allT_ATLonly_left);
aln_hld = aln_dta_sbj(strcmpi(aln_dta(:,2),'left'));
[ ~, ejk_nme, aln_nme ] = setxor( ejk_hld, aln_hld);
ejk_hld(ejk_nme)
aln_hld(aln_nme)

ejk_hld = cln_dta_sbj( grp.tle_post_allT_ATLonly_right);
aln_hld = aln_dta_sbj(strcmpi(aln_dta(:,2),'right'));
[ ~, ejk_nme, aln_nme ] = setxor( ejk_hld, aln_hld);
ejk_hld(ejk_nme)
aln_hld(aln_nme)

sbj_ind = find(strcmpi(cln_dta_sbj(:,1),'epd020'));
ismember(sbj_ind,intersect(tot_sbj,cog_sbj))
ismember(sbj_ind,intersect(intersect(tot_sbj,cog_sbj),str_sbj))
ismember(sbj_ind,intersect(intersect(intersect(tot_sbj,cog_sbj),str_sbj),srg_sbj))
ismember(sbj_ind,intersect(intersect(intersect(intersect(tot_sbj,cog_sbj),str_sbj),srg_sbj), lft_tle))
ismember(sbj_ind,intersect(intersect(intersect(intersect(tot_sbj,cog_sbj),str_sbj),srg_sbj), rgh_tle))

%% Compare wmparc-Entorhinal
aln_ind = find(strcmpi(aln_dta(:,2),'left'));
new_ind = grp.tle_post_3T_ATLonly_left;

[ ~, ovr_sbj_238, ovr_sbj_aln ] = intersect(cln_dta_sbj(new_ind), aln_dta_sbj(aln_ind));

[ aln_dta_sbj(aln_ind(ovr_sbj_aln)) aln_dta(aln_ind(ovr_sbj_aln), strcmpi(aln_dta_col,'fiber_FA-L_Unc')) ...
  cln_dta_sbj(new_ind(ovr_sbj_238)) fib_dta(new_ind(ovr_sbj_238), strcmpi(fib_dta_col,'L_Unc')) ]

