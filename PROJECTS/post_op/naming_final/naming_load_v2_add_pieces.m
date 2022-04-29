%% Get participants
% Load Scores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.red_fle = red_cap_fle;
[~ , ~ , ~ , sbj_cog, ~] = mmil_load_redcap(fcfg);

% Check who has post operative scores
sbj_nme_hld = cell(0);
for iT = cog_tst_inc
    sbj_nme_hld = [ sbj_nme_hld ; sbj_cog.sbj_nme(~isnan(sbj_cog.(cog_tst_nme{iT}))) ];
end
sbj_nme = unique(sbj_nme_hld);

clear sbj_cog

% Re-Load Redcap
fcfg = [];
fcfg.sbj_nme = sbj_nme;
fcfg.red_fle = red_cap_fle;
[sbj_dem , sbj_sze , sbj_scn , sbj_cog, sbj_emo, sbj_srg] = mmil_load_redcap(fcfg);

%% Raw post-op scores
ejk_chk_dir([prj_dir '/' prj_nme '/' 'Data' '/'])

cog_tst_nme = { 'bnt_raw_scr' 'ant_mem_raw_scr' 'bnt_raw_scr_pst' 'ant_mem_raw_scr_pst' 'bnt_raw_scr_pst_raw' 'ant_mem_raw_scr_pst_raw' };
cog_tst_col = { 'bnt_raw_scr' 'ant_mem_raw_scr' 'bnt_raw_scr_pst' 'ant_mem_raw_scr_pst' 'bnt_raw_scr_pst'     'ant_mem_raw_scr_pst' };

% Calculate Change Scores
fcfg = [];
fcfg.rci = 1;
[ pst_cog_dta , pst_cog_cat ] = ejk_post_cognitive(fcfg,sbj_cog,sbj_scn);

% Cognitive Save
cog_dta_out(:,1) = pst_cog_dta.sbj_nme;
for iC = 1:numel(cog_tst_nme)
    if iC==3 || iC==4
        cog_dta_out(:,iC+1) = num2cell(pst_cog_dta.(cog_tst_col{iC}));
    else
        cog_dta_out(:,iC+1) = num2cell(sbj_cog.(cog_tst_col{iC}));
    end
end
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive' '_' 'with_post_raw' '.csv'], [ 'SubjID' cog_tst_nme ; cog_dta_out ]);

% QC
cog_sbj_cln = { 'epd_ucsf010'                                                                              'epd_ucsf018'                                                                       'fc078' };
cog_roi_cln = { {'bnt_raw_scr_pst' 'ant_mem_raw_scr_pst' 'bnt_raw_scr_pst_raw' 'ant_mem_raw_scr_pst_raw'} {'bnt_raw_scr' 'ant_mem_raw_scr' 'bnt_raw_scr_pst_raw' 'ant_mem_raw_scr_pst_raw' } {'bnt_raw_scr' 'ant_mem_raw_scr' } };

cog_dta = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive' '_' 'with_post_raw' '.csv']);

fcfg = [];

fcfg.dta     = cog_dta;
fcfg.dta_col = 2:size(cog_dta,2);
fcfg.sbj_col = 1;

fcfg.sbj_nme = cog_sbj_cln;
fcfg.roi_nme = cog_roi_cln;

cln_dta = ejk_clean_roi(fcfg);

cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive' '_' 'with_post_raw' '_' 'QC' '.csv'], cln_dta)


%% Timing
sbj_scn


