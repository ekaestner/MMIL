
out_hld = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/TLE_Controls_pre_pre/ant_mem_raw_scr';

%
dta_lhs = load([ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/Data' '/' 'surf_wmparc_fa_lhs_sm' '313' '.mat']); %
dta_rhs = load([ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/Data' '/' 'surf_wmparc_fa_rhs_sm' '313' '.mat']);

%
load( [ '/home/ekaestne/PROJECTS/OUTPUT' '/' 'PostOperative/Naming' '/' 'groups.mat' ] );

%
cog_dta_nme = [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/Data' '/' 'Cognitive'                '_' 'QC' '.csv'];;
cog_dta = mmil_readtext(cog_dta_nme);
    cog_dta_col = ejk_fix_column_names(cog_dta(1,2:end));
    cog_dta_sbj = cog_dta(2:end,1);
    cog_dta     = cell2mat(cog_dta(2:end,2:end));

%
wmp_dta_nme = [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/Data' '/' 'wmparc_FA_wm_aparc_annot' '_' 'QC' '.csv'];
wmp_dta = mmil_readtext(wmp_dta_nme);
    wmp_dta_col = ejk_fix_column_names(wmp_dta(1,5:end));
    wmp_dta_sbj = wmp_dta(2:end,1);
    wmp_dta     = wmp_dta(2:end,5:end);

%
[ lhs_prc_loc, lhs_prc_lbl, ~]=fs_read_annotation( [ '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer' '/' 'lh' '.aparc.annot' ] );
[ rhs_prc_loc, rhs_prc_lbl, ~]=fs_read_annotation( [ '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer' '/' 'rh' '.aparc.annot' ] );
    
%%
deg_fre = sum(strcmpi(grp_var.(cfg.grp_nme{iG}),fst_nme)) - 2;

% Cluster threshold
srf_hld = fs_read_surf('/home/ekaestne/PROJECTS/EXTERNAL/home/rh.pial');
srf_chr = fs_calc_triarea(srf_hld);

ttt = fs_load_mgh('/home/ekaestne/PROJECTS/EXTERNAL/home/thickness-sphere-sm2819-lh.mgz');
bad_ind = find(ttt==0);
gdd_ind = 1:numel(ttt);
gdd_ind(bad_ind) = [];

fcfg = [];

fcfg.nverts = numel(ttt)-numel(bad_ind);
fcfg.area   = sum(srf_chr.vertex_area(gdd_ind));
fcfg.fwhm   = 1.25 * sqrt(cfg.smt_stp);
fcfg.df     = deg_fre;
fcfg.alpha  = cfg.pvl_cls;
fcfg.pval   = cfg.pvl_chs;

cls_thr = fs_calc_cluster_thresh( fcfg.nverts , ...
    fcfg.area , ...
    fcfg.fwhm , ...
    fcfg.df , ...
    fcfg.alpha , ...
    fcfg.pval );

pvl_nme = num2str(roundsd(cfg.pvl_chs,3));
pvl_nme = pvl_nme(3:end);
cls_nme = num2str(round(cls_thr));

% p-value Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvl_lhs = load([ out_hld '/' 'pvalues_lhs_dep_var_' fst_nme '.mat' ]); %  '_' 'sm' '313'
pvl_lhs = pvl_lhs.pvalues;
pvl_rhs = load([ out_hld '/' 'pvalues_rhs_dep_var_' fst_nme '.mat' ]); % '_' 'sm' '313'
pvl_rhs = pvl_rhs.pvalues;

% em-mean Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rvl_dff_lhs = load([ out_hld '/' 'rvalues_lhs_dep_var_' fst_nme '.mat' ]); %  '_' 'sm' '313'
rvl_dff_lhs = rvl_dff_lhs.rvalues;
rvl_dff_rhs = load([ out_hld '/' 'rvalues_rhs_dep_var_' fst_nme '.mat' ]); % '_' 'sm' '313'
rvl_dff_rhs = rvl_dff_rhs.rvalues;

% pvalue-cluster %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.pvl_thr = cfg.pvl_chs;
fcfg.cls_thr = cls_thr; % mm^2
fcfg.fsr_sbj = 'fsaverage';
fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc/';
fcfg.hms     = 'lh';
pvl_lhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_lhs );

fcfg.hms     = 'rh';
pvl_rhs_cls = ejk_surface_pvalue_cluster( fcfg , pvl_rhs );

%
lhs_bdd_ind = find(pvl_lhs_cls>fcfg.pvl_thr);
rhs_bdd_ind = find(pvl_rhs_cls>fcfg.pvl_thr);

rvl_dff_lhs_cls = rvl_dff_lhs;
rvl_dff_rhs_cls = rvl_dff_rhs;

rvl_dff_lhs_cls(lhs_bdd_ind) = 0;
rvl_dff_rhs_cls(rhs_bdd_ind) = 0;

%%
lhs_new_fus_roi = find(rvl_dff_lhs_cls>0);
rhs_new_fus_roi = find(rvl_dff_rhs_cls>0);

lhs_fus = find(lhs_prc_loc==8);
rhs_fus = find(rhs_prc_loc==8);

lhs_ovr_lap_fus = intersect(lhs_new_fus_roi,lhs_fus);
rhs_ovr_lap_fus = intersect(rhs_new_fus_roi,rhs_fus);

lhs_ovr_non_fus = setxor(lhs_fus,lhs_ovr_lap_fus);
rhs_ovr_non_fus = setxor(rhs_fus,rhs_ovr_lap_fus);


new_roi_hld = cell( size(dta_lhs.srf_dta_sbj,1), 3);

for iS = 1:size(dta_lhs.srf_dta_sbj,1)
    
    new_roi_hld{iS,1} = dta_lhs.srf_dta_sbj{iS};
    
    new_roi_hld{iS,2} = nanmean(dta_lhs.srf_dta(iS,lhs_new_fus_roi));
    new_roi_hld{iS,3} = nanmean(dta_rhs.srf_dta(iS,rhs_new_fus_roi));
    
    new_roi_hld{iS,4} = nanmean(dta_lhs.srf_dta(iS,lhs_fus));
    new_roi_hld{iS,5} = nanmean(dta_rhs.srf_dta(iS,rhs_fus));
    
    new_roi_hld{iS,6} = nanmean(dta_lhs.srf_dta(iS,lhs_ovr_lap_fus));
    new_roi_hld{iS,7} = nanmean(dta_rhs.srf_dta(iS,rhs_ovr_lap_fus));
    
    new_roi_hld{iS,8} = nanmean(dta_lhs.srf_dta(iS,lhs_ovr_non_fus));
    new_roi_hld{iS,9} = nanmean(dta_rhs.srf_dta(iS,rhs_ovr_non_fus));
    
end

%% PRE-ANT
cog_dta_hld = cog_dta(grp.tle_controls_pre_3T_allSurg_all, 2);
    cog_ind = find(~isnan(cog_dta_hld));
    
hld_old_lhs = cell2mat(wmp_dta(grp.tle_controls_pre_3T_allSurg_all,6));
    old_lhs_ind = find(~isnan(hld_old_lhs));
hld_old_rhs = cell2mat(wmp_dta(grp.tle_controls_pre_3T_allSurg_all,39));
    old_rhs_ind = find(~isnan(hld_old_rhs)); 
    
hld_new_tot_lhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,2));
    new_tot_lhs_ind = find(~isnan(hld_new_tot_lhs));
hld_new_tot_rhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,3));
    new_tot_rhs_ind = find(~isnan(hld_new_tot_rhs));

hld_new_fus_lhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,4));
    new_fus_lhs_ind = find(~isnan(hld_new_fus_lhs));
hld_new_fus_rhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,5));
    new_fus_rhs_ind = find(~isnan(hld_new_fus_rhs));

hld_new_ovr_lap_lhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,6));
    new_ovr_lap_lhs_ind = find(~isnan(hld_new_ovr_lap_lhs));
hld_new_ovr_lap_rhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,7));
    new_ovr_lap_rhs_ind = find(~isnan(hld_new_ovr_lap_rhs));

hld_new_ovr_non_lhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,8));
    new_ovr_non_lhs_ind = find(~isnan(hld_new_ovr_non_lhs));
hld_new_ovr_non_rhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,9));
    new_ovr_non_rhs_ind = find(~isnan(hld_new_ovr_non_rhs));
    
figure()

% OLD LHS - ANT
[ ~, pvl ] = corrcoef(cog_dta_hld( intersect(old_lhs_ind,cog_ind) ), hld_old_lhs( intersect(old_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr(cog_dta_hld( intersect(old_lhs_ind,cog_ind) ), hld_old_lhs( intersect(old_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,1)
scatter( cog_dta_hld( intersect(old_lhs_ind,cog_ind) ), hld_old_lhs( intersect(old_lhs_ind,cog_ind) ))
title(['OLD-LHS/ANT p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - TOTAL ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_tot_lhs_ind,cog_ind) ), hld_new_tot_lhs( intersect(new_tot_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_tot_lhs_ind,cog_ind) ), hld_new_tot_lhs( intersect(new_tot_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,5)
scatter( cog_dta_hld( intersect(new_tot_lhs_ind,cog_ind) ), hld_new_tot_lhs( intersect(new_tot_lhs_ind,cog_ind) ) )
title(['NEW-LHS/ANT TOTAL p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - FUS ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_fus_lhs_ind,cog_ind) ), hld_new_fus_lhs( intersect(new_fus_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_fus_lhs_ind,cog_ind) ), hld_new_fus_lhs( intersect(new_fus_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,6)
scatter( cog_dta_hld( intersect(new_fus_lhs_ind,cog_ind) ), hld_new_fus_lhs( intersect(new_fus_lhs_ind,cog_ind) ) )
title(['NEW-LHS/ANT FUS p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - OVR_LAP ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_lap_lhs_ind,cog_ind) ), hld_new_ovr_lap_lhs( intersect(new_ovr_lap_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_lap_lhs_ind,cog_ind) ), hld_new_ovr_lap_lhs( intersect(new_ovr_lap_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,7)
scatter( cog_dta_hld( intersect(new_ovr_lap_lhs_ind,cog_ind) ), hld_new_ovr_lap_lhs( intersect(new_ovr_lap_lhs_ind,cog_ind) ) )
title(['NEW-LHS/ANT Ovr-Lap p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - ovr_non ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_non_lhs_ind,cog_ind) ), hld_new_ovr_non_lhs( intersect(new_ovr_non_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_non_lhs_ind,cog_ind) ), hld_new_ovr_non_lhs( intersect(new_ovr_non_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,8)
scatter( cog_dta_hld( intersect(new_ovr_non_lhs_ind,cog_ind) ), hld_new_ovr_non_lhs( intersect(new_ovr_non_lhs_ind,cog_ind) ) )
title(['NEW-LHS/ANT Non-Ovr p=' num2str(roundsd(pvl(1,1),2))])

% OLD RHS - ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(old_rhs_ind,cog_ind) ), hld_old_rhs( intersect(old_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(old_rhs_ind,cog_ind) ), hld_old_rhs( intersect(old_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,9)
scatter( cog_dta_hld( intersect(old_rhs_ind,cog_ind) ), hld_old_rhs( intersect(old_rhs_ind,cog_ind) ))
title(['OLD-RHS/ANT p=' num2str(roundsd(pvl(1,1),2))])

% NEW RHS - ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_tot_rhs_ind,cog_ind) ), hld_new_tot_rhs( intersect(new_tot_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_tot_rhs_ind,cog_ind) ), hld_new_tot_rhs( intersect(new_tot_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,13)
scatter( cog_dta_hld( intersect(new_tot_rhs_ind,cog_ind) ), hld_new_tot_rhs( intersect(new_tot_rhs_ind,cog_ind) ) )
title(['NEW-rhs/ANT TOTAL p=' num2str(roundsd(pvl(1,1),2))])

% NEW rhs - FUS ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_fus_rhs_ind,cog_ind) ), hld_new_fus_rhs( intersect(new_fus_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_fus_rhs_ind,cog_ind) ), hld_new_fus_rhs( intersect(new_fus_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,14)
scatter( cog_dta_hld( intersect(new_fus_rhs_ind,cog_ind) ), hld_new_fus_rhs( intersect(new_fus_rhs_ind,cog_ind) ) )
title(['NEW-rhs/ANT FUS p=' num2str(roundsd(pvl(1,1),2))])

% NEW rhs - OVR_LAP ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_lap_rhs_ind,cog_ind) ), hld_new_ovr_lap_rhs( intersect(new_ovr_lap_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_lap_rhs_ind,cog_ind) ), hld_new_ovr_lap_rhs( intersect(new_ovr_lap_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,15)
scatter( cog_dta_hld( intersect(new_ovr_lap_rhs_ind,cog_ind) ), hld_new_ovr_lap_rhs( intersect(new_ovr_lap_rhs_ind,cog_ind) ) )
title(['NEW-rhs/ANT Ovr-Lap p=' num2str(roundsd(pvl(1,1),2))])

% NEW rhs - ovr_non ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_non_rhs_ind,cog_ind) ), hld_new_ovr_non_rhs( intersect(new_ovr_non_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_non_rhs_ind,cog_ind) ), hld_new_ovr_non_rhs( intersect(new_ovr_non_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,16)
scatter( cog_dta_hld( intersect(new_ovr_non_rhs_ind,cog_ind) ), hld_new_ovr_non_rhs( intersect(new_ovr_non_rhs_ind,cog_ind) ) )
title(['NEW-rhs/ANT Non-Ovr p=' num2str(roundsd(pvl(1,1),2))])

set(gcf,'Position',[0 0 1440 1080])
tightfig();
set(gcf,'Position',[0 0 1440 1080])
print('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/NewFusInvestgiation/ANT_Pre','-dpng')
close all

%% POST-ANT
cog_dta_hld = cog_dta(grp.tle_post_3T_ATLonly_left, 5);
    cog_ind = find(~isnan(cog_dta_hld));
  
hld_old_lhs = cell2mat(wmp_dta(grp.tle_post_3T_ATLonly_left,6));
    old_lhs_ind = find(~isnan(hld_old_lhs));
hld_old_rhs = cell2mat(wmp_dta(grp.tle_post_3T_ATLonly_left,39));
    old_rhs_ind = find(~isnan(hld_old_rhs)); 
    
hld_new_tot_lhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,2));
    new_tot_lhs_ind = find(~isnan(hld_new_tot_lhs));
hld_new_tot_rhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,3));
    new_tot_rhs_ind = find(~isnan(hld_new_tot_rhs));

hld_new_fus_lhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,4));
    new_fus_lhs_ind = find(~isnan(hld_new_fus_lhs));
hld_new_fus_rhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,5));
    new_fus_rhs_ind = find(~isnan(hld_new_fus_rhs));

hld_new_ovr_lap_lhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,6));
    new_ovr_lap_lhs_ind = find(~isnan(hld_new_ovr_lap_lhs));
hld_new_ovr_lap_rhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,7));
    new_ovr_lap_rhs_ind = find(~isnan(hld_new_ovr_lap_rhs));

hld_new_ovr_non_lhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,8));
    new_ovr_non_lhs_ind = find(~isnan(hld_new_ovr_non_lhs));
hld_new_ovr_non_rhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,9));
    new_ovr_non_rhs_ind = find(~isnan(hld_new_ovr_non_rhs));    
    
figure()

% OLD LHS - ANT
[ ~, pvl ] = corrcoef(cog_dta_hld( intersect(old_lhs_ind,cog_ind) ), hld_old_lhs( intersect(old_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr(cog_dta_hld( intersect(old_lhs_ind,cog_ind) ), hld_old_lhs( intersect(old_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,1)
scatter( cog_dta_hld( intersect(old_lhs_ind,cog_ind) ), hld_old_lhs( intersect(old_lhs_ind,cog_ind) ))
title(['OLD-LHS/ANT p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - TOTAL ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_tot_lhs_ind,cog_ind) ), hld_new_tot_lhs( intersect(new_tot_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_tot_lhs_ind,cog_ind) ), hld_new_tot_lhs( intersect(new_tot_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,5)
scatter( cog_dta_hld( intersect(new_tot_lhs_ind,cog_ind) ), hld_new_tot_lhs( intersect(new_tot_lhs_ind,cog_ind) ) )
title(['NEW-LHS/ANT TOTAL p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - FUS ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_fus_lhs_ind,cog_ind) ), hld_new_fus_lhs( intersect(new_fus_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_fus_lhs_ind,cog_ind) ), hld_new_fus_lhs( intersect(new_fus_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,6)
scatter( cog_dta_hld( intersect(new_fus_lhs_ind,cog_ind) ), hld_new_fus_lhs( intersect(new_fus_lhs_ind,cog_ind) ) )
title(['NEW-LHS/ANT FUS p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - OVR_LAP ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_lap_lhs_ind,cog_ind) ), hld_new_ovr_lap_lhs( intersect(new_ovr_lap_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_lap_lhs_ind,cog_ind) ), hld_new_ovr_lap_lhs( intersect(new_ovr_lap_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,7)
scatter( cog_dta_hld( intersect(new_ovr_lap_lhs_ind,cog_ind) ), hld_new_ovr_lap_lhs( intersect(new_ovr_lap_lhs_ind,cog_ind) ) )
title(['NEW-LHS/ANT Ovr-Lap p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - ovr_non ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_non_lhs_ind,cog_ind) ), hld_new_ovr_non_lhs( intersect(new_ovr_non_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_non_lhs_ind,cog_ind) ), hld_new_ovr_non_lhs( intersect(new_ovr_non_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,8)
scatter( cog_dta_hld( intersect(new_ovr_non_lhs_ind,cog_ind) ), hld_new_ovr_non_lhs( intersect(new_ovr_non_lhs_ind,cog_ind) ) )
title(['NEW-LHS/ANT Non-Ovr p=' num2str(roundsd(pvl(1,1),2))])

% OLD RHS - ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(old_rhs_ind,cog_ind) ), hld_old_rhs( intersect(old_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(old_rhs_ind,cog_ind) ), hld_old_rhs( intersect(old_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,9)
scatter( cog_dta_hld( intersect(old_rhs_ind,cog_ind) ), hld_old_rhs( intersect(old_rhs_ind,cog_ind) ))
title(['OLD-RHS/ANT p=' num2str(roundsd(pvl(1,1),2))])

% NEW RHS - ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_tot_rhs_ind,cog_ind) ), hld_new_tot_rhs( intersect(new_tot_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_tot_rhs_ind,cog_ind) ), hld_new_tot_rhs( intersect(new_tot_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,13)
scatter( cog_dta_hld( intersect(new_tot_rhs_ind,cog_ind) ), hld_new_tot_rhs( intersect(new_tot_rhs_ind,cog_ind) ) )
title(['NEW-rhs/ANT TOTAL p=' num2str(roundsd(pvl(1,1),2))])

% NEW rhs - FUS ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_fus_rhs_ind,cog_ind) ), hld_new_fus_rhs( intersect(new_fus_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_fus_rhs_ind,cog_ind) ), hld_new_fus_rhs( intersect(new_fus_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,14)
scatter( cog_dta_hld( intersect(new_fus_rhs_ind,cog_ind) ), hld_new_fus_rhs( intersect(new_fus_rhs_ind,cog_ind) ) )
title(['NEW-rhs/ANT FUS p=' num2str(roundsd(pvl(1,1),2))])

% NEW rhs - OVR_LAP ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_lap_rhs_ind,cog_ind) ), hld_new_ovr_lap_rhs( intersect(new_ovr_lap_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_lap_rhs_ind,cog_ind) ), hld_new_ovr_lap_rhs( intersect(new_ovr_lap_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,15)
scatter( cog_dta_hld( intersect(new_ovr_lap_rhs_ind,cog_ind) ), hld_new_ovr_lap_rhs( intersect(new_ovr_lap_rhs_ind,cog_ind) ) )
title(['NEW-rhs/ANT Ovr-Lap p=' num2str(roundsd(pvl(1,1),2))])

% NEW rhs - ovr_non ANT
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_non_rhs_ind,cog_ind) ), hld_new_ovr_non_rhs( intersect(new_ovr_non_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_non_rhs_ind,cog_ind) ), hld_new_ovr_non_rhs( intersect(new_ovr_non_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,16)
scatter( cog_dta_hld( intersect(new_ovr_non_rhs_ind,cog_ind) ), hld_new_ovr_non_rhs( intersect(new_ovr_non_rhs_ind,cog_ind) ) )
title(['NEW-rhs/ANT Non-Ovr p=' num2str(roundsd(pvl(1,1),2))])

set(gcf,'Position',[0 0 1440 1080])
tightfig();
set(gcf,'Position',[0 0 1440 1080])
print('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/NewFusInvestgiation/ANT_Post','-dpng')
close all

%% PRE-bnt
cog_dta_hld = cog_dta(grp.tle_controls_pre_3T_allSurg_all, 1);
    cog_ind = find(~isnan(cog_dta_hld));
    
hld_old_lhs = cell2mat(wmp_dta(grp.tle_controls_pre_3T_allSurg_all,6));
    old_lhs_ind = find(~isnan(hld_old_lhs));
hld_old_rhs = cell2mat(wmp_dta(grp.tle_controls_pre_3T_allSurg_all,39));
    old_rhs_ind = find(~isnan(hld_old_rhs)); 
    
hld_new_tot_lhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,2));
    new_tot_lhs_ind = find(~isnan(hld_new_tot_lhs));
hld_new_tot_rhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,3));
    new_tot_rhs_ind = find(~isnan(hld_new_tot_rhs));

hld_new_fus_lhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,4));
    new_fus_lhs_ind = find(~isnan(hld_new_fus_lhs));
hld_new_fus_rhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,5));
    new_fus_rhs_ind = find(~isnan(hld_new_fus_rhs));

hld_new_ovr_lap_lhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,6));
    new_ovr_lap_lhs_ind = find(~isnan(hld_new_ovr_lap_lhs));
hld_new_ovr_lap_rhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,7));
    new_ovr_lap_rhs_ind = find(~isnan(hld_new_ovr_lap_rhs));

hld_new_ovr_non_lhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,8));
    new_ovr_non_lhs_ind = find(~isnan(hld_new_ovr_non_lhs));
hld_new_ovr_non_rhs = cell2mat(new_roi_hld(grp.tle_controls_pre_3T_allSurg_all,9));
    new_ovr_non_rhs_ind = find(~isnan(hld_new_ovr_non_rhs));
    
figure()

% OLD LHS - bnt
[ ~, pvl ] = corrcoef(cog_dta_hld( intersect(old_lhs_ind,cog_ind) ), hld_old_lhs( intersect(old_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr(cog_dta_hld( intersect(old_lhs_ind,cog_ind) ), hld_old_lhs( intersect(old_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,1)
scatter( cog_dta_hld( intersect(old_lhs_ind,cog_ind) ), hld_old_lhs( intersect(old_lhs_ind,cog_ind) ))
title(['OLD-LHS/bnt p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - TOTAL bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_tot_lhs_ind,cog_ind) ), hld_new_tot_lhs( intersect(new_tot_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_tot_lhs_ind,cog_ind) ), hld_new_tot_lhs( intersect(new_tot_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,5)
scatter( cog_dta_hld( intersect(new_tot_lhs_ind,cog_ind) ), hld_new_tot_lhs( intersect(new_tot_lhs_ind,cog_ind) ) )
title(['NEW-LHS/bnt TOTAL p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - FUS bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_fus_lhs_ind,cog_ind) ), hld_new_fus_lhs( intersect(new_fus_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_fus_lhs_ind,cog_ind) ), hld_new_fus_lhs( intersect(new_fus_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,6)
scatter( cog_dta_hld( intersect(new_fus_lhs_ind,cog_ind) ), hld_new_fus_lhs( intersect(new_fus_lhs_ind,cog_ind) ) )
title(['NEW-LHS/bnt FUS p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - OVR_LAP bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_lap_lhs_ind,cog_ind) ), hld_new_ovr_lap_lhs( intersect(new_ovr_lap_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_lap_lhs_ind,cog_ind) ), hld_new_ovr_lap_lhs( intersect(new_ovr_lap_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,7)
scatter( cog_dta_hld( intersect(new_ovr_lap_lhs_ind,cog_ind) ), hld_new_ovr_lap_lhs( intersect(new_ovr_lap_lhs_ind,cog_ind) ) )
title(['NEW-LHS/bnt Ovr-Lap p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - ovr_non bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_non_lhs_ind,cog_ind) ), hld_new_ovr_non_lhs( intersect(new_ovr_non_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_non_lhs_ind,cog_ind) ), hld_new_ovr_non_lhs( intersect(new_ovr_non_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,8)
scatter( cog_dta_hld( intersect(new_ovr_non_lhs_ind,cog_ind) ), hld_new_ovr_non_lhs( intersect(new_ovr_non_lhs_ind,cog_ind) ) )
title(['NEW-LHS/bnt Non-Ovr p=' num2str(roundsd(pvl(1,1),2))])

% OLD RHS - bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(old_rhs_ind,cog_ind) ), hld_old_rhs( intersect(old_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(old_rhs_ind,cog_ind) ), hld_old_rhs( intersect(old_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,9)
scatter( cog_dta_hld( intersect(old_rhs_ind,cog_ind) ), hld_old_rhs( intersect(old_rhs_ind,cog_ind) ))
title(['OLD-RHS/bnt p=' num2str(roundsd(pvl(1,1),2))])

% NEW RHS - bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_tot_rhs_ind,cog_ind) ), hld_new_tot_rhs( intersect(new_tot_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_tot_rhs_ind,cog_ind) ), hld_new_tot_rhs( intersect(new_tot_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,13)
scatter( cog_dta_hld( intersect(new_tot_rhs_ind,cog_ind) ), hld_new_tot_rhs( intersect(new_tot_rhs_ind,cog_ind) ) )
title(['NEW-rhs/bnt TOTAL p=' num2str(roundsd(pvl(1,1),2))])

% NEW rhs - FUS bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_fus_rhs_ind,cog_ind) ), hld_new_fus_rhs( intersect(new_fus_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_fus_rhs_ind,cog_ind) ), hld_new_fus_rhs( intersect(new_fus_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,14)
scatter( cog_dta_hld( intersect(new_fus_rhs_ind,cog_ind) ), hld_new_fus_rhs( intersect(new_fus_rhs_ind,cog_ind) ) )
title(['NEW-rhs/bnt FUS p=' num2str(roundsd(pvl(1,1),2))])

% NEW rhs - OVR_LAP bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_lap_rhs_ind,cog_ind) ), hld_new_ovr_lap_rhs( intersect(new_ovr_lap_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_lap_rhs_ind,cog_ind) ), hld_new_ovr_lap_rhs( intersect(new_ovr_lap_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,15)
scatter( cog_dta_hld( intersect(new_ovr_lap_rhs_ind,cog_ind) ), hld_new_ovr_lap_rhs( intersect(new_ovr_lap_rhs_ind,cog_ind) ) )
title(['NEW-rhs/bnt Ovr-Lap p=' num2str(roundsd(pvl(1,1),2))])

% NEW rhs - ovr_non bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_non_rhs_ind,cog_ind) ), hld_new_ovr_non_rhs( intersect(new_ovr_non_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_non_rhs_ind,cog_ind) ), hld_new_ovr_non_rhs( intersect(new_ovr_non_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,16)
scatter( cog_dta_hld( intersect(new_ovr_non_rhs_ind,cog_ind) ), hld_new_ovr_non_rhs( intersect(new_ovr_non_rhs_ind,cog_ind) ) )
title(['NEW-rhs/bnt Non-Ovr p=' num2str(roundsd(pvl(1,1),2))])

set(gcf,'Position',[0 0 1440 1080])
tightfig();
set(gcf,'Position',[0 0 1440 1080])
print('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/NewFusInvestgiation/BNT_Pre','-dpng')
close all

%% POST-bnt
cog_dta_hld = cog_dta(grp.tle_post_3T_ATLonly_left, 4);
    cog_ind = find(~isnan(cog_dta_hld));
  
hld_old_lhs = cell2mat(wmp_dta(grp.tle_post_3T_ATLonly_left,6));
    old_lhs_ind = find(~isnan(hld_old_lhs));
hld_old_rhs = cell2mat(wmp_dta(grp.tle_post_3T_ATLonly_left,39));
    old_rhs_ind = find(~isnan(hld_old_rhs)); 
    
hld_new_tot_lhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,2));
    new_tot_lhs_ind = find(~isnan(hld_new_tot_lhs));
hld_new_tot_rhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,3));
    new_tot_rhs_ind = find(~isnan(hld_new_tot_rhs));

hld_new_fus_lhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,4));
    new_fus_lhs_ind = find(~isnan(hld_new_fus_lhs));
hld_new_fus_rhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,5));
    new_fus_rhs_ind = find(~isnan(hld_new_fus_rhs));

hld_new_ovr_lap_lhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,6));
    new_ovr_lap_lhs_ind = find(~isnan(hld_new_ovr_lap_lhs));
hld_new_ovr_lap_rhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,7));
    new_ovr_lap_rhs_ind = find(~isnan(hld_new_ovr_lap_rhs));

hld_new_ovr_non_lhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,8));
    new_ovr_non_lhs_ind = find(~isnan(hld_new_ovr_non_lhs));
hld_new_ovr_non_rhs = cell2mat(new_roi_hld(grp.tle_post_3T_ATLonly_left,9));
    new_ovr_non_rhs_ind = find(~isnan(hld_new_ovr_non_rhs));    
    
figure()

% OLD LHS - bnt
[ ~, pvl ] = corrcoef(cog_dta_hld( intersect(old_lhs_ind,cog_ind) ), hld_old_lhs( intersect(old_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr(cog_dta_hld( intersect(old_lhs_ind,cog_ind) ), hld_old_lhs( intersect(old_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,1)
scatter( cog_dta_hld( intersect(old_lhs_ind,cog_ind) ), hld_old_lhs( intersect(old_lhs_ind,cog_ind) ))
title(['OLD-LHS/bnt p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - TOTAL bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_tot_lhs_ind,cog_ind) ), hld_new_tot_lhs( intersect(new_tot_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_tot_lhs_ind,cog_ind) ), hld_new_tot_lhs( intersect(new_tot_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,5)
scatter( cog_dta_hld( intersect(new_tot_lhs_ind,cog_ind) ), hld_new_tot_lhs( intersect(new_tot_lhs_ind,cog_ind) ) )
title(['NEW-LHS/bnt TOTAL p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - FUS bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_fus_lhs_ind,cog_ind) ), hld_new_fus_lhs( intersect(new_fus_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_fus_lhs_ind,cog_ind) ), hld_new_fus_lhs( intersect(new_fus_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,6)
scatter( cog_dta_hld( intersect(new_fus_lhs_ind,cog_ind) ), hld_new_fus_lhs( intersect(new_fus_lhs_ind,cog_ind) ) )
title(['NEW-LHS/bnt FUS p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - OVR_LAP bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_lap_lhs_ind,cog_ind) ), hld_new_ovr_lap_lhs( intersect(new_ovr_lap_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_lap_lhs_ind,cog_ind) ), hld_new_ovr_lap_lhs( intersect(new_ovr_lap_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,7)
scatter( cog_dta_hld( intersect(new_ovr_lap_lhs_ind,cog_ind) ), hld_new_ovr_lap_lhs( intersect(new_ovr_lap_lhs_ind,cog_ind) ) )
title(['NEW-LHS/bnt Ovr-Lap p=' num2str(roundsd(pvl(1,1),2))])

% NEW LHS - ovr_non bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_non_lhs_ind,cog_ind) ), hld_new_ovr_non_lhs( intersect(new_ovr_non_lhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_non_lhs_ind,cog_ind) ), hld_new_ovr_non_lhs( intersect(new_ovr_non_lhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,8)
scatter( cog_dta_hld( intersect(new_ovr_non_lhs_ind,cog_ind) ), hld_new_ovr_non_lhs( intersect(new_ovr_non_lhs_ind,cog_ind) ) )
title(['NEW-LHS/bnt Non-Ovr p=' num2str(roundsd(pvl(1,1),2))])

% OLD RHS - bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(old_rhs_ind,cog_ind) ), hld_old_rhs( intersect(old_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(old_rhs_ind,cog_ind) ), hld_old_rhs( intersect(old_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,9)
scatter( cog_dta_hld( intersect(old_rhs_ind,cog_ind) ), hld_old_rhs( intersect(old_rhs_ind,cog_ind) ))
title(['OLD-RHS/bnt p=' num2str(roundsd(pvl(1,1),2))])

% NEW RHS - bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_tot_rhs_ind,cog_ind) ), hld_new_tot_rhs( intersect(new_tot_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_tot_rhs_ind,cog_ind) ), hld_new_tot_rhs( intersect(new_tot_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,13)
scatter( cog_dta_hld( intersect(new_tot_rhs_ind,cog_ind) ), hld_new_tot_rhs( intersect(new_tot_rhs_ind,cog_ind) ) )
title(['NEW-rhs/bnt TOTAL p=' num2str(roundsd(pvl(1,1),2))])

% NEW rhs - FUS bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_fus_rhs_ind,cog_ind) ), hld_new_fus_rhs( intersect(new_fus_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_fus_rhs_ind,cog_ind) ), hld_new_fus_rhs( intersect(new_fus_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,14)
scatter( cog_dta_hld( intersect(new_fus_rhs_ind,cog_ind) ), hld_new_fus_rhs( intersect(new_fus_rhs_ind,cog_ind) ) )
title(['NEW-rhs/bnt FUS p=' num2str(roundsd(pvl(1,1),2))])

% NEW rhs - OVR_LAP bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_lap_rhs_ind,cog_ind) ), hld_new_ovr_lap_rhs( intersect(new_ovr_lap_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_lap_rhs_ind,cog_ind) ), hld_new_ovr_lap_rhs( intersect(new_ovr_lap_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,15)
scatter( cog_dta_hld( intersect(new_ovr_lap_rhs_ind,cog_ind) ), hld_new_ovr_lap_rhs( intersect(new_ovr_lap_rhs_ind,cog_ind) ) )
title(['NEW-rhs/bnt Ovr-Lap p=' num2str(roundsd(pvl(1,1),2))])

% NEW rhs - ovr_non bnt
[ ~, pvl ] = corrcoef( cog_dta_hld( intersect(new_ovr_non_rhs_ind,cog_ind) ), hld_new_ovr_non_rhs( intersect(new_ovr_non_rhs_ind,cog_ind) ));
[ ~, pvl ] = corr( cog_dta_hld( intersect(new_ovr_non_rhs_ind,cog_ind) ), hld_new_ovr_non_rhs( intersect(new_ovr_non_rhs_ind,cog_ind) ), 'type', 'Spearman');
subplot(4,4,16)
scatter( cog_dta_hld( intersect(new_ovr_non_rhs_ind,cog_ind) ), hld_new_ovr_non_rhs( intersect(new_ovr_non_rhs_ind,cog_ind) ) )
title(['NEW-rhs/bnt Non-Ovr p=' num2str(roundsd(pvl(1,1),2))])

set(gcf,'Position',[0 0 1440 1080])
tightfig();
set(gcf,'Position',[0 0 1440 1080])
print('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/SurfaceCorrelation/NewFusInvestgiation/BNT_Post','-dpng')
close all





