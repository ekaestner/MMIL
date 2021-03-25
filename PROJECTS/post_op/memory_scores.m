lod_dta = mmil_readtext( [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena' '/' 'Memory_FA.csv' ] );
    lod_dta( cellfun(@isempty,lod_dta) ) = {NaN}; 
    
neu_bio_dta     = lod_dta(:,[1 13:36]);
neu_bio_sbj_nme = neu_bio_dta(2:end,1);
neu_bio_roi_nme = neu_bio_dta(1,2:end);
neu_bio_dta     = cell2mat(neu_bio_dta(2:end,2:end));

cog_scr_dta     = lod_dta(:,[1 5:12]);
cog_scr_sbj_nme = cog_scr_dta(2:end,1);
cog_scr_roi_nme = cog_scr_dta(1,2:end);
cog_scr_dta     = cell2mat(cog_scr_dta(2:end,2:end));

cov_dta     = lod_dta(:,1:4);
cov_sbj_nme = cov_dta(2:end,1);
cov_roi_nme = cov_dta(1,2:end);
cov_dta     = cov_dta(2:end,2:end);

sde_col = find(strcmpi(cov_roi_nme, 'SideOfSeizureFocus'));
lft_ind = find(strcmpi( cov_dta(:,sde_col), 'left'));
rgh_ind = find(strcmpi( cov_dta(:,sde_col), 'right'));

%%
% LEFT SIDE
fcfg = [];

fcfg.sbj_nme = neu_bio_sbj_nme(lft_ind,:);

fcfg.dta_one = cog_scr_dta(lft_ind,:);
fcfg.lbl_one = cog_scr_roi_nme;

fcfg.cor_typ = 'spearman';

fcfg.dta_two = neu_bio_dta(lft_ind,:);
fcfg.lbl_two = neu_bio_roi_nme;

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.15;

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/scores_by_wmparc_LEFT';

ejk_cross_cor( fcfg );

% RIGHT SIDE
fcfg = [];

fcfg.sbj_nme = neu_bio_sbj_nme;

fcfg.dta_one = neu_bio_dta;
fcfg.lbl_one = neu_bio_roi_nme;

fcfg.cor_typ = 'spearman';

fcfg.dta_two = cog_scr_dta;
fcfg.lbl_two = cog_scr_roi_nme;

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.15;

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/scores_by_wmparc_run';

ejk_cross_cor( fcfg );
