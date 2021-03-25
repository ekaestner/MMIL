prj_dir = '/home/ekaestne/PROJECTS/OUTPUT';
prj_nme = 'PostOperative/Naming';

clear out_put

nme_typ = { 'left_pre' 'right_pre' 'left_post' 'right_post'};
str_typ = 2; % 1: 1.5T & 3T; 2: 3T
    str_nme = { 'allT' '3T' };
srg_typ = 2; % 1: All; 2: ATL-only
    srg_nme = { 'allSurg' 'ATLonly' };

tst_pre = [1 2 3];
tst_pst = [4 5 6];

%%
cln_fle     = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv' ];

cog_fle     = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive.csv'];
cog_tst_nme = { 'bnt_nor_scr'     'ant_mem_raw_scr'     'cat_flu_nor_scr'     ...
                'bnt_nor_scr_pst' 'ant_mem_raw_scr_pst' 'cat_flu_nor_scr_pst' };

mri_mse     = { 'cort_thick_ctx' 'subcort_vol_ICV_cor' };
mri_roi_use = { [ 2 ]            [ 0 ] };

fmr_mse     = { 'alicia' };
fmr_roi_use = { [0]      };

dti_mse     = { 'wmparc_FA_wm' 'fiber_FA' };
dti_roi_use = { [ 2 ]          [ 0 ]      };

rsf_mse     = { 'var_ctx' 'var_vol' };
rsf_roi_use = { [ 2 ]     [ 0 ]     };

roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/';
roi_nme = { 'aparc.a2009s.annot' 'aparc.annot' };

%% Setup Groups
cln_dta = mmil_readtext(cln_fle);
cln_dta_col = cln_dta(1,2:end);
cln_dta_sbj = cln_dta(2:end,1);
cln_dta     = cln_dta(2:end,2:end);

cog_dta     = mmil_readtext(cog_fle);
cog_dta_col = cog_dta(1,2:end);
cog_dta_sbj = cog_dta(2:end,1);
cog_dta     = cog_dta(2:end,2:end);

% Choose subjects to use %%%%%%%%
% Find pre_patients
tot_pre = zeros(0);
for iT = 1:numel(tst_pre)
    tot_pre = [ tot_pre ; find(~isnan(cell2mat(cog_dta(:,tst_pre(iT))))) ];
end
tot_pre = unique(tot_pre);

% Find post_patients
tot_pst = zeros(0);
for iT = 1:numel(tst_pst)
    tot_pst = [ tot_pst ; find(~isnan(cell2mat(cog_dta(:,tst_pst(iT))))) ];
end
tot_pst = unique(tot_pst);

% Find side of epilepsy onset
lft_tle = find( strcmpi(cln_dta(:,2),'L') );
rgh_tle = find( strcmpi(cln_dta(:,2),'R') );

% Find Strength
if str_typ==1
    out_put = [prj_dir '/' prj_nme '/' 'CrossCorrelation' '_' str_nme{str_typ}];
elseif str_typ==2
    lft_tle = intersect( lft_tle, find(strcmpi( cellfun(@num2str,cln_dta(:,3),'uni',0), '3')) );
    rgh_tle = intersect( rgh_tle, find(strcmpi( cellfun(@num2str,cln_dta(:,3),'uni',0), '3')) );
    out_put = [prj_dir '/' prj_nme '/' 'CrossCorrelation' '_' str_nme{str_typ}];
end

% Find Surgery
if srg_typ==1
    lft_tle_srg = intersect( lft_tle, find(~cellfun(@isempty, cln_dta(:,10))) );
    rgh_tle_srg = intersect( rgh_tle, find(~cellfun(@isempty, cln_dta(:,10))) );    
    out_put = [ out_put '_' srg_nme{srg_typ} ];
elseif srg_typ==2
    lft_tle_srg = intersect( lft_tle, find(strcmpi( cln_dta(:,10), 'ATL')) );
    rgh_tle_srg = intersect( rgh_tle, find(strcmpi( cln_dta(:,10), 'ATL')) );
    out_put = [ out_put '_' srg_nme{srg_typ} ];
end

out_put = [ out_put '_' 'noQC'];

%% DTI
ts3_dti.(nme_typ{1}) = intersect(lft_tle, tot_pre) ; % 
ts3_dti.(nme_typ{2}) = intersect(rgh_tle, tot_pre) ; %
ts3_dti.(nme_typ{3}) = intersect(lft_tle_srg, tot_pst) ; %
ts3_dti.(nme_typ{4}) = intersect(rgh_tle_srg, tot_pst) ; %

cog_col.(nme_typ{1}) = tst_pre;
cog_col.(nme_typ{2}) = tst_pre;
cog_col.(nme_typ{3}) = [ tst_pre tst_pst ];
cog_col.(nme_typ{4}) = [ tst_pre tst_pst ];

for iG = 1:numel(nme_typ)
    for iN = 1:numel(dti_mse)
        for iR = 1:numel(dti_roi_use{iN})
            
            if ~(dti_roi_use{iN}(iR)==0)
                dti_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iN} '_' mmil_spec_char(roi_nme{dti_roi_use{iN}(iR)},{'.'}) '.csv'];
                dti_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iN} '_' mmil_spec_char(roi_nme{dti_roi_use{iN}(iR)},{'.'}) '_LI.csv'];
            else
                dti_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iN} '.csv'];
                dti_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iN} '_LI.csv'];
            end
            
            dti_dta = mmil_readtext(dti_dta_nme);
            dti_dta_col = ejk_fix_column_names(dti_dta(1,5:end));
            dti_dta_sbj = dti_dta(2:end,1);
            dti_dta     = dti_dta(2:end,5:end);
            dti_dta_lat = mmil_readtext(dti_dta_lat_nme);
            dti_dta_lat_col = dti_dta_lat(1,5:end);
            dti_dta_lat_sbj = dti_dta_lat(2:end,1);
            dti_dta_lat     = dti_dta_lat(2:end,5:end);
            
            % Raw
            fcfg = [];
            
            fcfg.sbj_nme = cln_dta_sbj( ts3_dti.(nme_typ{iG}), 1);
            
            fcfg.dta_one = cell2mat(dti_dta( ts3_dti.(nme_typ{iG}), :));
            fcfg.lbl_one = dti_dta_col;
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_two = cell2mat(cog_dta( ts3_dti.(nme_typ{iG}), cog_col.(nme_typ{iG})));
            fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.20;
            
            fcfg.out_dir = [ out_put '/' 'DTI' '/' dti_mse{iN} '/' 'Raw' '/' nme_typ{iG} '/'];
            
            ejk_cross_cor( fcfg );
            
            % LI
            fcfg = [];
            
            fcfg.sbj_nme = cln_dta_sbj( ts3_dti.(nme_typ{iG}), 1);
            
            fcfg.dta_one = cell2mat(dti_dta_lat( ts3_dti.(nme_typ{iG}), :));
            fcfg.lbl_one = dti_dta_lat_col;
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_two = cell2mat(cog_dta( ts3_dti.(nme_typ{iG}), cog_col.(nme_typ{iG})));
            fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.20;
            
            fcfg.out_dir = [ out_put '/' 'DTI' '/' dti_mse{iN} '/' 'LI' '/' nme_typ{iG} '/'];
            
            ejk_cross_cor( fcfg );
            
        end
    end
end

%% MRI
ts3_mri.(nme_typ{1}) = intersect(lft_tle, tot_pre) ; % 
ts3_mri.(nme_typ{2}) = intersect(rgh_tle, tot_pre) ; %
ts3_mri.(nme_typ{3}) = intersect(lft_tle_srg, tot_pst) ; %
ts3_mri.(nme_typ{4}) = intersect(rgh_tle_srg, tot_pst) ; %

cog_col.(nme_typ{1}) = tst_pre;
cog_col.(nme_typ{2}) = tst_pre;
cog_col.(nme_typ{3}) = [ tst_pre tst_pst ];
cog_col.(nme_typ{4}) = [ tst_pre tst_pst ];

for iG = 1:numel(nme_typ)
    for iN = 1:numel(mri_mse)
        for iR = 1:numel(mri_roi_use{iN})
            
            if ~(mri_roi_use{iN}(iR)==0)
                mri_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iN} '_' mmil_spec_char(roi_nme{mri_roi_use{iN}(iR)},{'.'}) '.csv'];
                mri_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iN} '_' mmil_spec_char(roi_nme{mri_roi_use{iN}(iR)},{'.'}) '_LI.csv'];
            else
                mri_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iN} '.csv'];
                mri_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iN} '_LI.csv'];
            end
            
            mri_dta = mmil_readtext(mri_dta_nme);
            mri_dta_col = ejk_fix_column_names(mri_dta(1,5:end));
            mri_dta_sbj = mri_dta(2:end,1);
            mri_dta     = mri_dta(2:end,5:end);
            mri_dta_lat = mmil_readtext(mri_dta_lat_nme);
            mri_dta_lat_col = mri_dta_lat(1,5:end);
            mri_dta_lat_sbj = mri_dta_lat(2:end,1);
            mri_dta_lat     = mri_dta_lat(2:end,5:end);
            
            % Raw
            fcfg = [];
            
            fcfg.sbj_nme = cln_dta_sbj( ts3_mri.(nme_typ{iG}), 1);
            
            fcfg.dta_one = cell2mat(mri_dta( ts3_mri.(nme_typ{iG}), :));
            fcfg.lbl_one = mri_dta_col;
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_two = cell2mat(cog_dta( ts3_mri.(nme_typ{iG}), cog_col.(nme_typ{iG})));
            fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.20;
            
            fcfg.out_dir = [ out_put '/' 'MRI' '/' mri_mse{iN} '/' 'Raw' '/' nme_typ{iG} '/'];
            
            ejk_cross_cor( fcfg );
            
            % LI
            fcfg = [];
            
            fcfg.sbj_nme = cln_dta_sbj( ts3_mri.(nme_typ{iG}), 1);
            
            fcfg.dta_one = cell2mat(mri_dta_lat( ts3_mri.(nme_typ{iG}), :));
            fcfg.lbl_one = mri_dta_lat_col;
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_two = cell2mat(cog_dta( ts3_mri.(nme_typ{iG}), cog_col.(nme_typ{iG})));
            fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.20;
            
            fcfg.out_dir = [ out_put '/' 'MRI' '/' mri_mse{iN} '/' 'LI' '/' nme_typ{iG} '/'];
            
            ejk_cross_cor( fcfg );
            
        end
    end
end

%% fMRI
ts3_fmr.(nme_typ{1}) = intersect(lft_tle, tot_pre) ; % 
ts3_fmr.(nme_typ{2}) = intersect(rgh_tle, tot_pre) ; %
ts3_fmr.(nme_typ{3}) = intersect(lft_tle_srg, tot_pst) ; %
ts3_fmr.(nme_typ{4}) = intersect(rgh_tle_srg, tot_pst) ; %

cog_col.(nme_typ{1}) = tst_pre;
cog_col.(nme_typ{2}) = tst_pre;
cog_col.(nme_typ{3}) = [ tst_pre tst_pst ];
cog_col.(nme_typ{4}) = [ tst_pre tst_pst ];

for iG = 1:numel(nme_typ)
    for iN = 1:numel(fmr_mse)
        for iR = 1:numel(fmr_roi_use{iN})
            
            if ~(fmr_roi_use{iN}(iR)==0)
                fmr_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' fmr_mse{iN} '_' mmil_spec_char(roi_nme{fmr_roi_use{iN}(iR)},{'.'}) '.csv'];
                fmr_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' fmr_mse{iN} '_' mmil_spec_char(roi_nme{fmr_roi_use{iN}(iR)},{'.'}) '_LI.csv'];
            else
                fmr_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' fmr_mse{iN} '.csv'];
                fmr_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' fmr_mse{iN} '_LI.csv'];
            end
            
            fmr_dta = mmil_readtext(fmr_dta_nme);
            fmr_dta_col = ejk_fix_column_names(fmr_dta(1,5:end));
            fmr_dta_sbj = fmr_dta(2:end,1);
            fmr_dta     = fmr_dta(2:end,5:end);
            fmr_dta_lat = mmil_readtext(fmr_dta_lat_nme);
            fmr_dta_lat_col = fmr_dta_lat(1,5:end);
            fmr_dta_lat_sbj = fmr_dta_lat(2:end,1);
            fmr_dta_lat     = fmr_dta_lat(2:end,5:end);
            
            % Raw
            fcfg = [];
            
            fcfg.sbj_nme = cln_dta_sbj( ts3_fmr.(nme_typ{iG}), 1);
            
            fcfg.dta_one = cell2mat(fmr_dta( ts3_fmr.(nme_typ{iG}), :));
            fcfg.lbl_one = fmr_dta_col;
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_two = cell2mat(cog_dta( ts3_fmr.(nme_typ{iG}), cog_col.(nme_typ{iG})));
            fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.20;
            
            fcfg.out_dir = [ out_put '/' 'fMRI' '/' fmr_mse{iN} '/' 'Raw' '/' nme_typ{iG} '/'];
            
            ejk_cross_cor( fcfg );
            
            % LI
            fcfg = [];
            
            fcfg.sbj_nme = cln_dta_sbj( ts3_fmr.(nme_typ{iG}), 1);
            
            fcfg.dta_one = cell2mat(fmr_dta_lat( ts3_fmr.(nme_typ{iG}), :));
            fcfg.lbl_one = fmr_dta_lat_col;
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_two = cell2mat(cog_dta( ts3_fmr.(nme_typ{iG}), cog_col.(nme_typ{iG})));
            fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.20;
            
            fcfg.out_dir = [ out_put '/' 'fMRI' '/' fmr_mse{iN} '/' 'LI' '/' nme_typ{iG} '/'];
            
            ejk_cross_cor( fcfg );
            
        end
    end
end

%% rsfMRI
ts3_rsf.(nme_typ{1}) = intersect(lft_tle, tot_pre) ; % 
ts3_rsf.(nme_typ{2}) = intersect(rgh_tle, tot_pre) ; %
ts3_rsf.(nme_typ{3}) = intersect(lft_tle_srg, tot_pst) ; %
ts3_rsf.(nme_typ{4}) = intersect(rgh_tle_srg, tot_pst) ; %

cog_col.(nme_typ{1}) = tst_pre;
cog_col.(nme_typ{2}) = tst_pre;
cog_col.(nme_typ{3}) = [ tst_pre tst_pst ];
cog_col.(nme_typ{4}) = [ tst_pre tst_pst ];

for iG = 1:numel(nme_typ)
    for iN = 1:numel(rsf_mse)
        for iR = 1:numel(rsf_roi_use{iN})
            
            if ~(rsf_roi_use{iN}(iR)==0)
                rsf_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iN} '_' mmil_spec_char(roi_nme{rsf_roi_use{iN}(iR)},{'.'}) '.csv'];
                rsf_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iN} '_' mmil_spec_char(roi_nme{rsf_roi_use{iN}(iR)},{'.'}) '_LI.csv'];
            else
                rsf_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iN} '.csv'];
                rsf_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iN} '_LI.csv'];
            end
            
            rsf_dta = mmil_readtext(rsf_dta_nme);
            rsf_dta_col = ejk_fix_column_names(rsf_dta(1,5:end));
            rsf_dta_sbj = rsf_dta(2:end,1);
            rsf_dta     = rsf_dta(2:end,5:end);
            rsf_dta_lat = mmil_readtext(rsf_dta_lat_nme);
            rsf_dta_lat_col = rsf_dta_lat(1,5:end);
            rsf_dta_lat_sbj = rsf_dta_lat(2:end,1);
            rsf_dta_lat     = rsf_dta_lat(2:end,5:end);
            
            % Raw
            fcfg = [];
            
            fcfg.sbj_nme = cln_dta_sbj( ts3_rsf.(nme_typ{iG}), 1);
            
            fcfg.dta_one = cell2mat(rsf_dta( ts3_rsf.(nme_typ{iG}), :));
            fcfg.lbl_one = rsf_dta_col;
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_two = cell2mat(cog_dta( ts3_rsf.(nme_typ{iG}), cog_col.(nme_typ{iG})));
            fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.20;
            
            fcfg.out_dir = [ out_put '/' 'rsfMRI' '/' rsf_mse{iN} '/' 'Raw' '/' nme_typ{iG} '/'];
            
            ejk_cross_cor( fcfg );
            
            % LI
            fcfg = [];
            
            fcfg.sbj_nme = cln_dta_sbj( ts3_rsf.(nme_typ{iG}), 1);
            
            fcfg.dta_one = cell2mat(rsf_dta_lat( ts3_rsf.(nme_typ{iG}), :));
            fcfg.lbl_one = rsf_dta_lat_col;
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_two = cell2mat(cog_dta( ts3_rsf.(nme_typ{iG}), cog_col.(nme_typ{iG})));
            fcfg.lbl_two = cog_dta_col(cog_col.(nme_typ{iG}));
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.20;
            
            fcfg.out_dir = [ out_put '/' 'rsfMRI' '/' rsf_mse{iN} '/' 'LI' '/' nme_typ{iG} '/'];
            
            ejk_cross_cor( fcfg );
            
        end
    end
end

