clear; clc;

prj_dir = '/home/ekaestne/PROJECTS/OUTPUT';
prj_nme = 'PostOperative/Naming_final_sample';

%% Load Data
% Cognitive Data
fcfg = [];
fcfg.dta_loc = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Collaborations/Alena/nme_aml/CFScores.csv';
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

emp_hld = cell(1,numel(cog_dta_col));
for iC = 1:numel(cog_dta_col)
    if isnumeric(cog_dta{1,iC})
        emp_hld{1,iC} = NaN;
    else
        emp_hld{1,iC} = '';
    end
end

% Wmparc
wmp_dta     = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'wmparc_FA_wm' '_' 'aparc_annot' '_QC.csv']);
wmp_dta_col = ejk_fix_column_names(wmp_dta(1,5:end));
wmp_dta_sbj = wmp_dta(2:end,1);
wmp_dta     = wmp_dta(2:end,5:end);

% Fibers
trc_dta     = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'fiber_FA' '_QC.csv']);
trc_dta_col = ejk_fix_column_names(trc_dta(1,5:end));
trc_dta_sbj = trc_dta(2:end,1);
trc_dta     = trc_dta(2:end,5:end);

%% Setup Data for correlations
% Reorder Alena data
new_cog_dta     = cell(numel(wmp_dta_sbj), size(cog_dta,2));

for iS = 1:numel(wmp_dta_sbj)
    sbj_ind = find(strcmpi(cog_dta_sbj,wmp_dta_sbj{iS}));
    if isempty(sbj_ind)
        new_cog_dta(iS,:) = emp_hld;
    else
        new_cog_dta(iS,:) = cog_dta(sbj_ind,:);
    end    
end

% Create grps
grp.tle_lft = find( strcmpi(new_cog_dta(:,1),'Temporal Lobe Epilepsy Surg') & strcmpi(new_cog_dta(:,2),'left') );
grp.tle_rgh = find( strcmpi(new_cog_dta(:,1),'Temporal Lobe Epilepsy Surg') & strcmpi(new_cog_dta(:,2),'right') );
grp.con     = find( strcmpi(new_cog_dta(:,1),'Normal Control') );
grp.ttl     = sort([ grp.tle_lft ; grp.tle_rgh ; grp.con ]);

grp_nme = fieldnames(grp);

% Create Cognitive Data
cog_ind = ismember( cog_dta_col, { 'Pre_Category Fluency Total' 'Pre_Animals' 'Pre_FirstNames' } );
cog_cor_dta     = new_cog_dta(:,cog_ind);
cog_cor_dta_col = cog_dta_col(1,cog_ind);
    cog_cor_dta_col{1} = 'Pre_Category_Fluency_Total';

% Create post-operative Cognitive Data
cog_ind = ismember( cog_dta_col, { 'Post_Category Fluency Total' 'Post_Animals' 'Post_FirstNames' } );
cog_cor_dta_pst     = new_cog_dta(:,cog_ind);
    cog_cor_dta_pst = cog_cor_dta_pst - cog_cor_dta;
    cog_cor_dta_pst(cellfun(@isempty,cog_cor_dta_pst)) = {NaN};
cog_cor_dta_pst_col = cog_dta_col(1,cog_ind);
    cog_cor_dta_pst_col{1} = 'Post_Category_Fluency_Total';
    
% Create Neurobio Data
neu_cor_dta     = [ trc_dta     wmp_dta ];
neu_cor_dta_col = [ trc_dta_col wmp_dta_col ];

%% Run correlations
out_put = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Collaborations/Alena/nme_aml/correlations'; ejk_chk_dir(out_put);
for iG = 1:numel(grp_nme)
    
    fcfg = [];
    
    fcfg.sbj_nme = wmp_dta_sbj( grp.(grp_nme{iG}), 1);
    
    fcfg.dta_one = cell2mat(neu_cor_dta( grp.(grp_nme{iG}), :));
    fcfg.lbl_one = neu_cor_dta_col;
    
    fcfg.cor_typ = 'spearman';
    
    fcfg.dta_two = cell2mat(cog_cor_dta( grp.(grp_nme{iG}), :));
    fcfg.lbl_two = cog_cor_dta_col;
    
    fcfg.pvl_cut = 0.05;
    fcfg.pvl_lib = 0.10;
    
    fcfg.out_dir = [ out_put '/' grp_nme{iG} '/'];
    
    ejk_cross_cor( fcfg );
    
end

for iG = 1:numel(grp_nme)
    
    fcfg = [];
    
    fcfg.sbj_nme = wmp_dta_sbj( grp.(grp_nme{iG}), 1);
    
    fcfg.dta_one = cell2mat(neu_cor_dta( grp.(grp_nme{iG}), :));
    fcfg.lbl_one = neu_cor_dta_col;
    
    fcfg.cor_typ = 'spearman';
    
    fcfg.dta_two = cell2mat(cog_cor_dta_pst( grp.(grp_nme{iG}), :));
    fcfg.lbl_two = cog_cor_dta_pst_col;
    
    fcfg.pvl_cut = 0.05;
    fcfg.pvl_lib = 0.10;
    
    fcfg.out_dir = [ out_put '/' grp_nme{iG} '_post' '/'];
    
    ejk_cross_cor( fcfg );
    
end










