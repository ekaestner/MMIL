prj_dir = '/home/ekaestne/PROJECTS/OUTPUT';
prj_nme = 'PostOperative/Language';

run_typ = 3;

%%
cln_fle     = [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Language/Data' '/' 'Clinical.csv' ];

cog_fle     = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Language/Data/Cognitive_QC.csv';
cog_tst_nme = { 'bnt_nor_scr'     'ant_mem_raw_scr'     'cat_flu_nor_scr'     'ltr_tot_nor_scr' ...
    'bnt_nor_scr_pst' 'ant_mem_raw_scr_pst' 'cat_flu_nor_scr_pst' 'ltr_tot_nor_scr_pst' };

mri_mse     = { 'cort_thick_ctx' 'subcort_vol_ICV_cor' };
mri_roi_use = { [ 2 ]            [ 0 ] };

dti_mse     = { 'wmparc_FA_wm' 'fiber_FA' };
dti_roi_use = { [ 2 ]          [ 0 ]      };

rsf_mse     = { 'var_ctx' 'var_vol' };
rsf_roi_use = { [ 2 ]     [ 0 ]     };

roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/';
roi_nme = { 'aparc.a2009s.annot' 'aparc.annot' };

inc_sbj = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Language/Data/Inclusion.csv');
mri_inc = find( cell2mat(inc_sbj(2:end, strcmpi(inc_sbj(1,:),'MRI_Include'))));
dti_inc = find( cell2mat(inc_sbj(2:end, strcmpi(inc_sbj(1,:),'DTI_Include'))));
fmr_inc = find( cell2mat(inc_sbj(2:end, strcmpi(inc_sbj(1,:),'fMRI_Include'))));
rsf_inc = find( cell2mat(inc_sbj(2:end, strcmpi(inc_sbj(1,:),'rsfMRI_Include'))));

%% Setup Groups
cln_dta = mmil_readtext(cln_fle);
cln_dta_col = cln_dta(1,2:end);
cln_dta_sbj = cln_dta(2:end,1);
cln_dta     = cln_dta(2:end,2:end);

cog_dta     = mmil_readtext(cog_fle);
cog_dta_col = cog_dta(1,2:end);
cog_dta_sbj = cog_dta(2:end,1);
cog_dta     = cog_dta(2:end,2:end);

% Choose subjects to use
if run_typ==1
    out_put = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Language/CrossCorrelation_3T/';
    ts3_lft = find( strcmpi(cln_dta(:,2),'L') & strcmpi( cellfun(@num2str,cln_dta(:,3),'uni',0), '3') );
    ts3_rgh = find( strcmpi(cln_dta(:,2),'R') & strcmpi( cellfun(@num2str,cln_dta(:,3),'uni',0), '3') );
elseif run_typ==2
    out_put = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Language/CrossCorrelation_3T_ATL/';
    ts3_lft = find( strcmpi(cln_dta(:,2),'L') & strcmpi( cellfun(@num2str,cln_dta(:,3),'uni',0), '3') & strcmpi( cln_dta(:,10), 'ATL') );
    ts3_rgh = find( strcmpi(cln_dta(:,2),'R') & strcmpi( cellfun(@num2str,cln_dta(:,3),'uni',0), '3') & strcmpi( cln_dta(:,10), 'ATL') );
elseif run_typ==3
    out_put = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Language/CrossCorrelation_3T_ATL_Dominant/';
    ts3_lft = find( strcmpi(cln_dta(:,2),'L') & strcmpi( cellfun(@num2str,cln_dta(:,3),'uni',0), '3') & strcmpi( cln_dta(:,10), 'ATL') &  strcmpi( cln_dta(:,15), 'dominant') );
    ts3_rgh = find( strcmpi(cln_dta(:,2),'R') & strcmpi( cellfun(@num2str,cln_dta(:,3),'uni',0), '3') & strcmpi( cln_dta(:,10), 'ATL') &  strcmpi( cln_dta(:,15), 'dominant') );
end

%% DTI
ts3_lft_dti = intersect(ts3_lft, dti_inc);
ts3_rgh_dti = intersect(ts3_rgh, dti_inc);

for iN = 1:numel(dti_mse)
    for iR = 1:numel(dti_roi_use{iN})
        
        if ~(dti_roi_use{iN}(iR)==0)
            dti_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iN} '_' mmil_spec_char(roi_nme{dti_roi_use{iN}(iR)},{'.'}) '_QC.csv'];
            %             dti_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iN} '_' mmil_spec_char(roi_nme{dti_roi_use{iN}(iR)},{'.'}) '_LI_QC.csv'];
        else
            dti_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iN} '_QC.csv'];
            %             dti_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' dti_mse{iN} '_LI_QC.csv'];
        end
        
        dti_dta = mmil_readtext(dti_dta_nme);
        dti_dta_col = ejk_fix_column_names(dti_dta(1,5:end));
        dti_dta_sbj = dti_dta(2:end,1);
        dti_dta     = dti_dta(2:end,5:end);
        %         dti_dta_lat = mmil_readtext(dti_dta_lat_nme);
        %             dti_dta_lat_col = dti_dta_lat(1,5:end);
        %             dti_dta_lat_sbj = dti_dta_lat(2:end,1);
        %             dti_dta_lat     = dti_dta_lat(2:end,5:end);
        
        % Raw
        fcfg = [];
        
        fcfg.sbj_nme = cln_dta_sbj( ts3_lft_dti, 1);
        
        fcfg.dta_one = cell2mat(dti_dta( ts3_lft_dti, :));
        fcfg.lbl_one = dti_dta_col;
        
        fcfg.cor_typ = 'spearman';
        
        fcfg.dta_two = cell2mat(cog_dta( ts3_lft_dti, :));
        fcfg.lbl_two = cog_dta_col;
        
        fcfg.pvl_cut = 0.05;
        fcfg.pvl_lib = 0.20;
        
        fcfg.out_dir = [ out_put '/' 'DTI' '/' dti_mse{iN} '/' 'Raw' '/'];
        
        ejk_cross_cor( fcfg );
        
        % LI
        %         fcfg = [];
        %
        %         fcfg.sbj_nme = cln_dta_sbj( ts3_lft_dti, 1);
        %
        %         fcfg.dta_one = cell2mat(dti_dta_lat( ts3_lft_dti, :));
        %         fcfg.lbl_one = dti_dta_lat_col;
        %
        %         fcfg.cor_typ = 'spearman';
        %
        %         fcfg.dta_two = cell2mat(cog_dta( ts3_lft_dti, :));
        %         fcfg.lbl_two = cog_dta_col;
        %
        %         fcfg.pvl_cut = 0.05;
        %         fcfg.pvl_lib = 0.20;
        %
        %         fcfg.out_dir = [ out_put '/' 'DTI' '/' dti_mse{iN} '/' 'LI' '/'];
        %
        %         ejk_cross_cor( fcfg );
        
    end
end

%% MRI
ts3_lft_mri = intersect(ts3_lft, mri_inc);
ts3_rgh_mri = intersect(ts3_rgh, mri_inc);

for iN = 1:numel(mri_mse)
    for iR = 1:numel(mri_roi_use{iN})
        
        if ~(mri_roi_use{iN}(iR)==0)
            mri_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iN} '_' mmil_spec_char(roi_nme{mri_roi_use{iN}(iR)},{'.'}) '_QC.csv'];
            %             mri_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iN} '_' mmil_spec_char(roi_nme{mri_roi_use{iN}(iR)},{'.'}) '_LI.csv'];
        else
            mri_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iN} '_QC.csv'];
            %             mri_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' mri_mse{iN} '_LI.csv'];
        end
        
        mri_dta = mmil_readtext(mri_dta_nme);
        mri_dta_col = ejk_fix_column_names(mri_dta(1,5:end));
        mri_dta_sbj = mri_dta(2:end,1);
        mri_dta     = mri_dta(2:end,5:end);
        %         mri_dta_lat = mmil_readtext(mri_dta_lat_nme);
        %             mri_dta_lat_col = ejk_fix_column_names(mri_dta_lat(1,5:end));
        %             mri_dta_lat_sbj = mri_dta_lat(2:end,1);
        %             mri_dta_lat     = mri_dta_lat(2:end,5:end);
        
        % Raw
        fcfg = [];
        
        fcfg.sbj_nme = cln_dta_sbj( ts3_lft_mri, 1);
        
        fcfg.dta_one = cell2mat(mri_dta( ts3_lft_mri, :));
        fcfg.lbl_one = mri_dta_col;
        
        fcfg.cor_typ = 'spearman';
        
        fcfg.dta_two = cell2mat(cog_dta( ts3_lft_mri, :));
        fcfg.lbl_two = cog_dta_col;
        
        fcfg.pvl_cut = 0.05;
        fcfg.pvl_lib = 0.20;
        
        fcfg.out_dir = [ out_put '/' 'MRI' '/' mri_mse{iN} '/' 'Raw' '/'];
        
        ejk_cross_cor( fcfg );
        
        % LI
        %         fcfg = [];
        %
        %         fcfg.sbj_nme = cln_dta_sbj( ts3_lft_mri, 1);
        %
        %         fcfg.dta_one = cell2mat(mri_dta_lat( ts3_lft_mri, :));
        %         fcfg.lbl_one = mri_dta_lat_col;
        %
        %         fcfg.cor_typ = 'spearman';
        %
        %         fcfg.dta_two = cell2mat(cog_dta( ts3_lft_mri, :));
        %         fcfg.lbl_two = cog_dta_col;
        %
        %         fcfg.pvl_cut = 0.05;
        %         fcfg.pvl_lib = 0.20;
        %
        %         fcfg.out_dir = [ out_put '/' 'MRI' '/' mri_mse{iN} '/' 'LI' '/'];
        %
        %         ejk_cross_cor( fcfg );
        
    end
end

%% rsfMRI
ts3_lft_fmr = intersect(ts3_lft, rsf_inc);
ts3_rgh_fmr = intersect(ts3_rgh, rsf_inc);

for iN = 1:numel(rsf_mse)
    for iR = 1:numel(rsf_roi_use{iN})
        
        if ~(rsf_roi_use{iN}(iR)==0)
            rsf_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iN} '_' mmil_spec_char(roi_nme{rsf_roi_use{iN}(iR)},{'.'}) '_QC.csv'];
            %             rsf_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iN} '_' mmil_spec_char(roi_nme{rsf_roi_use{iN}(iR)},{'.'}) '_LI.csv'];
        else
            rsf_dta_nme     = [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iN} '_QC.csv'];
            %             rsf_dta_lat_nme = [prj_dir '/' prj_nme '/' 'Data' '/' rsf_mse{iN} '_LI.csv'];
        end
        
        rsf_dta = mmil_readtext(rsf_dta_nme);
        rsf_dta_col = ejk_fix_column_names(rsf_dta(1,5:end));
        rsf_dta_sbj = rsf_dta(2:end,1);
        rsf_dta     = rsf_dta(2:end,5:end);
        %         if iN~=3
        %         rsf_dta_lat = mmil_readtext(rsf_dta_lat_nme);
        %             rsf_dta_lat_col = ejk_fix_column_names(rsf_dta_lat(1,5:end));
        %             rsf_dta_lat_sbj = rsf_dta_lat(2:end,1);
        %             rsf_dta_lat     = rsf_dta_lat(2:end,5:end);
        %         end
        
        % Raw
        fcfg = [];
        
        fcfg.sbj_nme = cln_dta_sbj( ts3_lft_fmr, 1);
        
        fcfg.dta_one = cell2mat(rsf_dta( ts3_lft_fmr, :));
        fcfg.lbl_one = rsf_dta_col;
        
        fcfg.cor_typ = 'spearman';
        
        fcfg.dta_two = cell2mat(cog_dta( ts3_lft_fmr, :));
        fcfg.lbl_two = cog_dta_col;
        
        fcfg.pvl_cut = 0.05;
        fcfg.pvl_lib = 0.20;
        
        fcfg.out_dir = [ out_put '/' 'rsfMRI' '/' rsf_mse{iN} '/' 'Raw' '/'];
        
        ejk_cross_cor( fcfg );
        
        % LI
        %         if iN ~= 3
        %         fcfg = [];
        %
        %         fcfg.sbj_nme = cln_dta_sbj( ts3_lft_fmr, 1);
        %
        %         fcfg.dta_one = cell2mat(rsf_dta_lat( ts3_lft_fmr, :));
        %         fcfg.lbl_one = rsf_dta_lat_col;
        %
        %         fcfg.cor_typ = 'spearman';
        %
        %         fcfg.dta_two = cell2mat(cog_dta( ts3_lft_fmr, :));
        %         fcfg.lbl_two = cog_dta_col;
        %
        %         fcfg.pvl_cut = 0.05;
        %         fcfg.pvl_lib = 0.20;
        %
        %         fcfg.out_dir = [ out_put '/' 'rsfMRI' '/' rsf_mse{iN} '/' 'LI' '/'];
        %
        %         ejk_cross_cor( fcfg );
    end    
end

