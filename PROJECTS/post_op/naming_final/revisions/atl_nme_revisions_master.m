out_put = [ prj_dir '/' prj_nme '/' 'Revisions'];

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

%% Get WTAR, Language Psych, Duration, Test intervals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
red_cap_fle = '/home/ekaestne/PROJECTS/DATA/csv/redcap/Redcap_2021_11_17.csv';

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

% Re-Load Redcap
fcfg = [];
fcfg.sbj_nme = cog_dta_sbj;
fcfg.red_fle = red_cap_fle;
[sbj_dem , sbj_sze , sbj_scn , sbj_cog, sbj_emo, sbj_srg] = mmil_load_redcap(fcfg);

% Calculate Change Scores
fcfg = [];
fcfg.rci = 1;
[ pst_cog_dta , pst_cog_cat ] = ejk_post_cognitive(fcfg,sbj_cog,sbj_scn);

% New output data
rvs_dta_nme = { 'SubjID'     ...
                'Duration' ...
                'WTAR_pre' 'LF_pre' 'LF_post_rci' 'CF_pre' 'CF_post_rci' ...
                'PreInterval' 'PostInterval' };
rvs_dta_out = [ sbj_dem.sbj_nme        ...
                num2cell(sbj_sze.sbj_sze_dur)           ...
                num2cell(sbj_cog.wtr_int_nor_scr) num2cell(sbj_cog.ltr_tot_raw_scr) num2cell(pst_cog_dta.ltr_tot_nor_scr_pst) num2cell(sbj_cog.cat_flu_raw_scr) num2cell(pst_cog_dta.cat_flu_nor_scr_pst) ...
                num2cell(sbj_cog.neu_psy_tst_dte_gap) num2cell(sbj_cog.neu_psy_tst_dte_pst_gap) ];
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'Revision_Data' '.csv'], [ rvs_dta_nme ; rvs_dta_out ]);

%% Investigate Engel Scores
aln_dta = mmil_readtext( [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/AllData' '/' 'Memory_FA_Updated_Raw.csv' ] );
aln_dta = aln_dta(2:end,[1 16]);

red_cap_fle = '/home/ekaestne/PROJECTS/DATA/csv/redcap/Redcap_2021_11_17.csv';

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

% Re-Load Redcap
fcfg = [];
fcfg.sbj_nme = cog_dta_sbj;
fcfg.red_fle = red_cap_fle;
[sbj_dem , sbj_sze , sbj_scn , sbj_cog, sbj_emo, sbj_srg] = mmil_load_redcap(fcfg);

%
lft_eng = [ sbj_srg.sbj_nme( grp.tle_post_3T_ATLonly_left, 1) sbj_srg.eng_out( grp.tle_post_3T_ATLonly_left, :) ];
for iS = 1:size(lft_eng,1)
    ind = find(strcmpi(aln_dta(:,1), lft_eng{iS,1} ));
    if ~isempty(ind)
        lft_eng{iS,3} = aln_dta{ind,1};
        lft_eng{iS,4} = aln_dta{ind,2};
    else
        lft_eng{iS,3} = '';
        lft_eng{iS,4} = '';
    end
end

rgh_eng = [ sbj_srg.sbj_nme( grp.tle_post_3T_ATLonly_right, 1) sbj_srg.eng_out( grp.tle_post_3T_ATLonly_right, :) ]; 
for iS = 1:size(rgh_eng,1)
    ind = find(strcmpi(aln_dta(:,1), rgh_eng{iS,1} ));
    if ~isempty(ind)
        rgh_eng{iS,3} = aln_dta{ind,1};
        rgh_eng{iS,4} = aln_dta{ind,2};
    else
        rgh_eng{iS,3} = '';
        rgh_eng{iS,4} = '';
    end
end

% {'epd048' }{0×0 char}{'epd048' }{'II+'} - no simple engel, just extended engel
% {'epd_ucsf024'}{0×0 char}{'epd_ucsf024'}{'I'  } - no year 1 follow-up
% {'epd_ucsf026'}{0×0 char}{'epd_ucsf026'}{'I'  } - no year 1 follow-up
% {'epd_ucsf029'}{0×0 char}{'epd_ucsf029'}{'I'  } - no year 1 follow-up
% {'epd_ucsf031'}{0×0 char}{'epd_ucsf031'}{'I'  } - no year 1 follow-up
% {'epd_ucsf032'}{0×0 char}{'epd_ucsf032'}{'II+'} - no year 1 follow-up
% {'epd_ucsf037'}{0×0 char}{'epd_ucsf037'}{'II+'} - no year 1 follow-up
% {'epd_ucsf039'}{0×0 char}{'epd_ucsf039'}{'II+'} - no year 1 follow-up
% {'epd_ucsf042'}{0×0 char}{'epd_ucsf042'}{'I'  } - no year 1 follow-up

%% UCSD/UCSF overlap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_dir = [ out_put '/' 'Harmony' '/' ];
ejk_chk_dir(out_dir);

trc_dta     = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'fiber_FA' '_QC.csv']);
trc_dta_col = ejk_fix_column_names(trc_dta(1,5:end));
trc_dta_sbj = trc_dta(2:end,1);
trc_dta     = trc_dta(2:end,5:end);

wmp_dta     = mmil_readtext([prj_dir '/' prj_nme '/' 'Data' '/' 'wmparc_FA_wm' '_' 'aparc_annot' '_QC.csv']);
wmp_dta_col = ejk_fix_column_names(wmp_dta(1,5:end));
wmp_dta_sbj = wmp_dta(2:end,1);
wmp_dta     = wmp_dta(2:end,5:end);

[ ~, trc_dta_use ] = intersect( trc_dta_col, { 'xL_IFO' 'xR_IFO' 'xL_ILF' 'xR_ILF' });
[ ~, wmp_dta_use ] = intersect( wmp_dta_col, { 'xlh_fusiform' 'xrh_fusiform' });

% All TLE %%%%%%
san_frn_ind = string_find( trc_dta_sbj( grp.tle_pre_3T_allSurg_all ), 'ucsf');
san_dgo_ind = string_find( trc_dta_sbj( grp.tle_pre_3T_allSurg_all ), 'epd');
    san_dgo_ind = setxor( san_dgo_ind, san_frn_ind );
    
pvl_cnt = 1;
for iP = 1:numel(trc_dta_use)
    min_bin = min(cell2mat(trc_dta( grp.tle_pre_3T_allSurg_all, trc_dta_use(1))));
    max_bin = max(cell2mat(trc_dta( grp.tle_pre_3T_allSurg_all, trc_dta_use(1))));
    bin_edg = linspace(min_bin,max_bin,20);
    
    figure(); hold on;
    san_frn_dta = cell2mat(trc_dta( grp.tle_pre_3T_allSurg_all(san_frn_ind), trc_dta_use(iP)));
    san_dgo_dta = cell2mat(trc_dta( grp.tle_pre_3T_allSurg_all(san_dgo_ind), trc_dta_use(iP)));
    histogram( san_frn_dta, 'BinEdges',bin_edg);
    histogram( san_dgo_dta, 'BinEdges',bin_edg);
    max_num = get(gca,'ylim');
    line([nanmean(san_frn_dta) nanmean(san_frn_dta)],[0 max_num(2)+1],'LineWidth',4,'Color',rgb('blue'));
    line([nanmean(san_dgo_dta) nanmean(san_dgo_dta)],[0 max_num(2)],'LineWidth',4,'Color',rgb('orange'));
    set(gca,'xlim',[.30 .65]);
    [~,pvl] = ttest2(san_frn_dta,san_dgo_dta);
    pvl_hld(pvl_cnt) = pvl; pvl_cnt = pvl_cnt + 1;
%     title( [ trc_dta_col{trc_dta_use(iP)} '  ' num2str(roundsd(pvl,2)) ] );
    print(gcf,[ out_dir '/' trc_dta_col{trc_dta_use(iP)} '.png'],'-dpng');
    close all
end

for iP = 1:numel(wmp_dta_use)
    min_bin = min(cell2mat(wmp_dta( grp.tle_pre_3T_allSurg_all, wmp_dta_use(1))));
    max_bin = max(cell2mat(wmp_dta( grp.tle_pre_3T_allSurg_all, wmp_dta_use(1))));
    bin_edg = linspace(min_bin,max_bin,20);
    
    figure(); hold on;
    san_frn_dta = cell2mat(wmp_dta( grp.tle_pre_3T_allSurg_all(san_frn_ind), wmp_dta_use(iP)));
    san_dgo_dta = cell2mat(wmp_dta( grp.tle_pre_3T_allSurg_all(san_dgo_ind), wmp_dta_use(iP)));
    histogram( san_frn_dta, 'BinEdges',bin_edg);
    histogram( san_dgo_dta, 'BinEdges',bin_edg);
    max_num = get(gca,'ylim');
    line([nanmean(san_frn_dta) nanmean(san_frn_dta)],[0 max_num(2)+1],'LineWidth',4,'Color',rgb('blue'));
    line([nanmean(san_dgo_dta) nanmean(san_dgo_dta)],[0 max_num(2)],'LineWidth',4,'Color',rgb('orange'));
    set(gca,'xlim',[.15 .50]);
    [~,pvl] = ttest2(san_frn_dta,san_dgo_dta);
    pvl_hld(pvl_cnt) = pvl; pvl_cnt = pvl_cnt + 1;
%     title( [ trc_dta_col{trc_dta_use(iP)} '  ' num2str(roundsd(pvl,2)) ] );
    print(gcf,[ out_dir '/' wmp_dta_col{wmp_dta_use(iP)} '.png'],'-dpng');
    close all
end
cell2csv([out_dir '/' 'pvl_hld.csv'],[ num2cell(pvl_hld) median(pvl_hld) ])

%% Remove Trendlines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atl_nme_no_trend_line_figure

%% Additional correlations (LF, CF, Duration) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_dir = [ out_put '/' 'Correlations' '/' ];
ejk_chk_dir(out_dir);

% Load %%%%%%%%%%%
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Revision_Data_engel.csv'];
fcfg.dta_col = 2;
[ rev_dta, rev_dta_sbj, rev_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'subcort_vol_ICV_cor_QC.csv'];
fcfg.dta_col = 2;
[ hip_dta, hip_dta_sbj, hip_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'fiber_FA_QC.csv'];
fcfg.dta_col = 2;
[ trc_dta, trc_dta_sbj, trc_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'wmparc_FA_wm_aparc_annot_QC.csv'];
fcfg.dta_col = 2;
[ wmp_dta, wmp_dta_sbj, wmp_dta_col] = ejk_dta_frm( fcfg );

% Create Data 2 %%%%%%%%%%%
dta_two = cell2mat([ cog_dta(:,strcmpi(cog_dta_col,'bnt_raw_scr')) ...
                     cog_dta(:,strcmpi(cog_dta_col,'ant_mem_raw_scr')) ...
                     rev_dta(:,strcmpi(rev_dta_col,'LF_pre')) ...
                     rev_dta(:,strcmpi(rev_dta_col,'CF_pre')) ...
                     cog_dta(:,strcmpi(cog_dta_col,'bnt_raw_scr_pst')) ...
                     cog_dta(:,strcmpi(cog_dta_col,'ant_mem_raw_scr_pst')) ...
                     rev_dta(:,strcmpi(rev_dta_col,'LF_post_rci')) ...
                     rev_dta(:,strcmpi(rev_dta_col,'CF_post_rci')) ]);

dta_two_col = { 'bnt_raw_scr' ...
                'ant_mem_raw_scr' ...
                'LF_pre' ...
                'CF_pre' ...
                'bnt_raw_scr_pst' ...
                'ant_mem_raw_scr_pst' ...
                'LF_post_rci' ...
                'CF_post_rci' };

% Create Data 1 %%%%%%%%%%%
dta_one = cell2mat([ cln_dta(:,strcmpi(cln_dta_col,'Educ')) ...
                     cln_dta(:,strcmpi(cln_dta_col,'AgeAtSurgery')) ...
                     cln_dta(:,strcmpi(cln_dta_col,'AgeOfSeizureOnset')) ...
                     cln_dta(:,strcmpi(cln_dta_col,'NumAEDs')) ...
                     rev_dta(:,strcmpi(rev_dta_col,'Duration')) ...
                     rev_dta(:,strcmpi(rev_dta_col,'WTAR_pre')) ...
                     hip_dta(:,strcmpi(hip_dta_col,'Left_Hippocampus')) ...
                     hip_dta(:,strcmpi(hip_dta_col,'Right_Hippocampus')) ...
                     trc_dta(:,strcmpi(trc_dta_col,'L_ILF')) ...
                     trc_dta(:,strcmpi(trc_dta_col,'R_ILF')) ...
                     trc_dta(:,strcmpi(trc_dta_col,'L_IFO')) ...
                     trc_dta(:,strcmpi(trc_dta_col,'R_IFO')) ...
                     trc_dta(:,strcmpi(trc_dta_col,'L_tSLF')) ...
                     trc_dta(:,strcmpi(trc_dta_col,'R_tSLF')) ...
                     trc_dta(:,strcmpi(trc_dta_col,'L_Unc')) ...
                     trc_dta(:,strcmpi(trc_dta_col,'R_Unc')) ...
                     wmp_dta(:,strcmpi(wmp_dta_col,'lh_fusiform')) ...
                     wmp_dta(:,strcmpi(wmp_dta_col,'rh_fusiform')) ]);

dta_one_col = { 'Educ' ...
                'AgeAtSurgery' ...
                'AgeOfSeizureOnset' ...
                'NumAEDs' ...
                'Duration' ...
                'WTAR_pre' ...
                'Left_Hippocampus' ...
                'Right_Hippocampus' ...
                'L_ILF' ...
                'R_ILF' ...
                'L_IFO' ...
                'R_IFO' ...
                'L_tSLF' ...
                'R_tSLF' ...
                'L_Unc' ...
                'R_Unc' ...
                'lh_fusiform' ...
                'rh_fusiform' };
        
% Correlations for PRE %%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = cln_dta_sbj( grp.tle_controls_pre_3T_allSurg_all, 1);

fcfg.dta_one = dta_one( grp.tle_controls_pre_3T_allSurg_all, :);
fcfg.lbl_one = dta_one_col;

fcfg.cor_typ = 'spearman';

fcfg.dta_two = dta_two( grp.tle_controls_pre_3T_allSurg_all, 1:4);
fcfg.lbl_two = dta_two_col;

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_dir '/' 'PreSurgical' '/' ];

ejk_cross_cor( fcfg );

% Correlations for Post - LEFT %%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = cln_dta_sbj( grp.tle_post_3T_ATLonly_left, 1);

fcfg.dta_one = dta_one( grp.tle_post_3T_ATLonly_left, :);
fcfg.lbl_one = dta_one_col;

fcfg.cor_typ = 'spearman';

fcfg.dta_two = dta_two( grp.tle_post_3T_ATLonly_left, 5:8 );
fcfg.lbl_two = dta_two_col;

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_dir '/' 'PostSurgical_Left' '/' ];

ejk_cross_cor( fcfg );

% Correlations for Post - Right %%%%%%%%%%%
fcfg = [];

fcfg.sbj_nme = cln_dta_sbj( grp.tle_post_3T_ATLonly_right, 1);

fcfg.dta_one = dta_one( grp.tle_post_3T_ATLonly_right, :);
fcfg.lbl_one = dta_one_col;

fcfg.cor_typ = 'spearman';

fcfg.dta_two = dta_two( grp.tle_post_3T_ATLonly_right, 5:8 );
fcfg.lbl_two = dta_two_col;

fcfg.pvl_cut = 0.05;
fcfg.pvl_lib = 0.20;

fcfg.out_dir = [ out_dir '/' 'PostSurgical_Right' '/' ];

ejk_cross_cor( fcfg );

%% Additional Group Comparisons (Engel, Duration) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_dir = [ out_put '/' 'Group' '/' ];
ejk_chk_dir(out_dir);

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Revision_Data_engel.csv'];
fcfg.dta_col = 2;
[ rev_dta, rev_dta_sbj, rev_dta_col] = ejk_dta_frm( fcfg );

% t-test (Duration, WTAR) %%%%%%%%%%%
[~, use_dta_col ] = intersect( rev_dta_col, { 'Duration' 'WTAR_pre' });

fcfg = [];
fcfg.grp     = grp;
fcfg.grp_inc = {{'tle_post_3T_ATLonly_left' 'tle_post_3T_ATLonly_right'}};
fcfg.grp_nme = {{'LTLE'                    'RTLE' }};
fcfg.dta = cell2mat(rev_dta(:,use_dta_col));
fcfg.sbj = rev_dta_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );

% Test
fcfg = [];

fcfg.sbj_nme = grp_sbj{1};

fcfg.dta     = grp_dta{1};
fcfg.dta_nme = rev_dta_col(:,use_dta_col);

fcfg.grp     = grp_typ{1};
fcfg.grp_nme = {'tle_pre_3T_ATLonly_left_right'};

fcfg.out_dir = [ out_dir '/' 'ttest2' '/'];

ejk_ttest2_independent( fcfg );

% ANOVA (WTAR) %%%%%%%%%%%
[~, use_dta_col ] = intersect( rev_dta_col, { 'WTAR_pre' });

fcfg = [];
fcfg.grp     = grp;
fcfg.grp_inc = {{'controls_pre_3T_allSurg_all' 'tle_pre_3T_allSurg_left' 'tle_pre_3T_allSurg_right'}};
fcfg.grp_nme = {{'HC'                          'LTLE'                    'RTLE' }};
fcfg.dta = cell2mat(rev_dta(:,use_dta_col));
fcfg.sbj = rev_dta_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );

% Test
cfg.sbj_nme = grp_sbj{1};

cfg.dta     = grp_dta{1};
cfg.dta_nme = rev_dta_col(use_dta_col);

cfg.grp     = grp_typ{1};
cfg.grp_nme = {'tle_con_pre'};

cfg.out_dir = [ out_dir '/' 'ANOVA' '/'];

ejk_1way_anova( cfg )
    
% Fishers (Engel) %%%%%%%%%%%
[~, use_dta_col ] = intersect( rev_dta_col, { 'EngelOutcome' });

fcfg = [];
fcfg.grp     = grp;
fcfg.grp_inc = {{'tle_pre_3T_ATLonly_left' 'tle_pre_3T_ATLonly_right'}};
fcfg.grp_nme = {{'LTLE'                    'RTLE' }};
fcfg.dta = rev_dta(:,use_dta_col);
fcfg.sbj = rev_dta_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );

grp_dta{1}( cellfun(@isempty,grp_dta{1})) = {NaN};

% Test
fcfg = [];
            
fcfg.sbj = grp_sbj{1};

fcfg.dta_one = grp_dta{1};
fcfg.lbl_one = cln_dta_col(use_dta_col);

fcfg.dta_two = repmat(grp_typ{1},1,numel(use_dta_col));
fcfg.lbl_two = strcat( 'group_', cln_dta_col(use_dta_col));

fcfg.out_dir = [ out_dir '/' 'fishers' '/'];

ejk_fisher_test( fcfg );

%% Additionall Demo/Clinical Tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_dir = [ out_put '/' 'Tables' '/' ];
ejk_chk_dir(out_dir);

% Post
dta_inp{1} = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Revision_Data_engel.csv']);
dta_inp{1}( strcmpi(dta_inp{1}(:,10),'II') ,10) = {'II+'};
dta_inp{1}( strcmpi(dta_inp{1}(:,10),'III') ,10) = {'II+'};
dta_inp{1}( strcmpi(dta_inp{1}(:,10),'IV') ,10) = {'II+'};
dta_inp{2} = mmil_readtext([ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/Revisions/Group/fishers' '/' 'output_table.csv']);
dta_inp{3} = mmil_readtext([ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/Revisions/Group/ttest2/tle.pre.3T.ATLonly.left.right' '/' 'output_table.csv']);


fcfg = [];
fcfg.tbl(1,:) = { 'count,1,tle_post_3T_ATLonly_left,EngelOutcome,I/II+' ...
                  'count,1,tle_post_3T_ATLonly_right,EngelOutcome,I/II+' ...
                  'copy,2,SurgicalSide,report'};
fcfg.tbl(2,:) = { 'mean/std,1,tle_post_3T_ATLonly_left,Duration' ...
                  'mean/std,1,tle_post_3T_ATLonly_right,Duration' ...
                  'copy,3,Duration,report'};
fcfg.tbl(3,:) = { 'mean/std,1,tle_post_3T_ATLonly_left,WTAR_pre' ...
                  'mean/std,1,tle_post_3T_ATLonly_right,WTAR_pre' ...
                  'copy,3,WTAR.pre,report'};
fcfg.tbl(4,:) = { 'mean/std,1,tle_post_3T_ATLonly_left,PreInterval' ...
                  'mean/std,1,tle_post_3T_ATLonly_right,PreInterval' ...
                  'mean/std,1,tle_post_3T_ATLonly_all,PreInterval'};
fcfg.tbl(5,:) = { 'mean/std,1,tle_post_3T_ATLonly_left,PostInterval' ...
                  'mean/std,1,tle_post_3T_ATLonly_right,PostInterval' ...
                  'mean/std,1,tle_post_3T_ATLonly_all,PostInterval'};
fcfg.dta = dta_inp;
fcfg.grp = grp;
tbl_out = ejk_create_table( fcfg );

col_lbl = { '' 'L-TLE' 'R-TLE' 'Test -or- All-TLE' };
row_lbl = { 'Engel Outcome (I/II+)' 'Duration (years)' 'WTAR (pre)' 'Pre NeuroPsych -to- Imaging' 'Post NeuroPsych -to- Surgery' }';
out_tbl = [ col_lbl ; row_lbl tbl_out ];

cell2csv( [ out_dir '/' 'RevisionsTableOutput_post.csv' ], out_tbl)

% Pre
dta_inp{1} = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Revision_Data_engel.csv']);
dta_inp{1}( strcmpi(dta_inp{1}(:,10),'II') ,10) = {'II+'};
dta_inp{1}( strcmpi(dta_inp{1}(:,10),'III') ,10) = {'II+'};
dta_inp{1}( strcmpi(dta_inp{1}(:,10),'IV') ,10) = {'II+'};
dta_inp{2} = mmil_readtext([ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/Revisions/Group/ANOVA/tle.con.pre' '/' 'output_table.csv']);

fcfg = [];
fcfg.tbl(1,:) = { 'mean/std,1,controls_pre_3T_allSurg_all,WTAR_pre' ...
                  'mean/std,1,tle_post_3T_ATLonly_left,WTAR_pre' ...
                  'mean/std,1,tle_post_3T_ATLonly_right,WTAR_pre' ...
                  'copy,2,WTAR.pre,report'};

fcfg.dta = dta_inp;
fcfg.grp = grp;
tbl_out = ejk_create_table( fcfg );

col_lbl = { '' 'HC' 'L-TLE' 'R-TLE' 'Test -or- All-TLE' };
row_lbl = { 'WTAR (pre)' }';
out_tbl = [ col_lbl ; row_lbl tbl_out ];

cell2csv( [ out_dir '/' 'RevisionsTableOutput_pre.csv' ], out_tbl)

%% Additional Tract Correlation Table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear out_tbl fdr_use

out_dir = [ out_put '/' 'AdditionalCorrelations' '/' ];
ejk_chk_dir(out_dir);

dta_loc = { [ prj_dir '/' prj_nme '/' 'SpecificCor/DTI/fiber_FA/Raw/tle_controls_pre_3T_allSurg_all' '/' ] ...
            [ prj_dir '/' prj_nme '/' 'SpecificCor/DTI/fiber_FA/Raw/tle_post_3T_ATLonly_left' '/' ] ...
            [ prj_dir '/' prj_nme '/' 'SpecificCor/DTI/fiber_FA/Raw/tle_post_3T_ATLonly_right' '/' ] };

tbl_nme = { 'xL.UNC' 1 ; ...
            'xR.UNC' 1 ; ...
            'xL.tSLF' 1 ; ...
            'xR.tSLF' 1 ; 
            'xL.UNC' 2 ; ...
            'xR.UNC' 2 ; ...
            'xL.tSLF' 2 ; ...
            'xR.tSLF' 2 ; 
            'xL.UNC' 3 ; ...
            'xR.UNC' 3 ; ...
            'xL.tSLF' 3 ; ...
            'xR.tSLF' 3 };
        
col_lbl = { '' 'Pre-operative BNT' 'Pre-operative ANT' };
row_lbl = { 'Pre L-UNC' 2 ;...
            'Pre R-UNC' 2 ;...
            'Pre L-ARC' 2 ;...
            'Pre R-ARC' 2 ; ...
            'Post-LTLE L-UNC' 2 ;...
            'Post-LTLE R-UNC' 2 ;...
            'Post-LTLE L-ARC' 2 ;...
            'Post-LTLE R-ARC' 2
            'Post-RTLE L-UNC' 2 ;...
            'Post-RTLE R-UNC' 2 ;...
            'Post-RTLE L-ARC' 2 ;...
            'Post-RTLE R-ARC' 2 };

% Out Table - OUTPUT (r-value)
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_rvalues.csv' ] ); end

% %%%
fcfg = [];

for iR = 1:size(tbl_nme,1)
  if iR<5
  fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr'] ... 
                     ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',ant.mem.raw.scr'] };
  else
      fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr.pst'] ... 
                         ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',ant.mem.raw.scr.pst'] };
  end
end

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_out = ejk_create_table( fcfg );

% Supp TABLE 1 - Significance (p-value)
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_pvalues.csv' ] ); end

% %%%
fcfg = [];

for iR = 1:size(tbl_nme,1)
    if iR<5
  fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr'] ... 
                    ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',ant.mem.raw.scr'] };
    else
        fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr.pst'] ... 
                    ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',ant.mem.raw.scr.pst'] };
    end
end

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_pvl = ejk_create_table( fcfg );

% Supp TABLE 2 - Stat Out (summary)
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_stats.csv' ] ); end

% %%%
fcfg = [];

for iR = 1:size(tbl_nme,1)
if iR<5
    fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr'] ... 
                    ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',ant.mem.raw.scr'] };
else
    fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr.pst'] ... 
                    ['copy,' num2str(tbl_nme{iR,2}) ',' num2str(tbl_nme{iR,1}) ',ant.mem.raw.scr.pst'] };
end
end

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_stt = ejk_create_table( fcfg );

% Supp TABLE 3 - FDR (summary)
fdr_grp_hld = cell2mat(row_lbl(:,2));
fdr_typ_hld = unique(fdr_grp_hld);

for iF = 1:numel(fdr_typ_hld)
    fdr_grp_mbm = fdr_grp_hld==fdr_typ_hld(iF);
    fdr_hld = cell2mat(tbl_pvl(fdr_grp_mbm,:));
    fdr_use{iF} = FDR(fdr_hld(:),.05);
    if isempty(fdr_use{iF}); fdr_use{iF} = 0; end
end

tbl_fdr=tbl_pvl;
for iR = 1:size(tbl_out,1)
    fdr_ind = fdr_typ_hld==row_lbl{iR,2};
    for iC = 1:size(tbl_out,2)
       if tbl_fdr{iR,iC}>fdr_use{fdr_ind}
           tbl_fdr{iR,iC}=NaN;
       end       
    end
end

% Format table
for iR = 1:size(tbl_out,1)
    for iC = 1:size(tbl_out,2)
        out_tbl{iR,iC} = num2str(roundsd(tbl_out{iR,iC},2));
        
        if tbl_pvl{iR,iC} < .01
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '**'];
        elseif tbl_pvl{iR,iC} < .05
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '*'];
        end
        
        if ~isnan(tbl_fdr{iR,iC})
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '#'];
        end
        
    end
end

out_tbl = [ col_lbl ; row_lbl(:,1) out_tbl ];
out_stt = [ col_lbl ; row_lbl(:,1) tbl_stt ];

% Save Table
cell2csv( [ out_dir '/' 'AdditionalTracts.csv' ], out_tbl)
cell2csv( [ out_dir '/' 'AdditionalTracts_stats.csv' ], out_stt)

%% Additional Test Correlation Table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear out_tbl fdr_use

out_dir = [ out_put '/' 'AdditionalCorrelations' '/' ];
ejk_chk_dir(out_dir);

dta_loc = { [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/Revisions/Correlations/PreSurgical' ] ...
            ['/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/Revisions/Correlations/PostSurgical' ] };

tbl_nme = { 'Left.Hippocampus' 1 ;...
            'Right.Hippocampus' 1 ; ...
            'L.ILF' 1 ; ...
            'R.ILF' 1 ; ...
            'L.IFO' 1 ; ...
            'R.IFO' 1 ; ...
            'lh.fusiform' 1 ; ...
            'rh.fusiform' 1 };
        
col_lbl = { '' 'Pre-operative CF' 'Pre-operative LF' 'LTLE Post-operative CF' 'LTLE Post-operative LF' };
row_lbl = { 'Left.Hippocampus' 2 ;...
            'Right.Hippocampus' 2; ...
            'Left ILF' 2; ...
            'Right ILF' 2; ...
            'Left IFOF' 2; ...
            'Right IFOF' 2; ...
            'Left Fusiform' 2; ...
            'Right Fusiform' 2 };

% Out Table - OUTPUT (r-value)
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_rvalues.csv' ] ); end

% %%%
fcfg = [];

for iR = 1:size(tbl_nme,1)

  fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2})   ',' num2str(tbl_nme{iR,1}) ',CF.pre'] ... 
                     ['copy,' num2str(tbl_nme{iR,2})   ',' num2str(tbl_nme{iR,1}) ',LF.pre'] ...
                     ['copy,' num2str(tbl_nme{iR,2}+1) ',' num2str(tbl_nme{iR,1}) ',CF.pre'] ... 
                     ['copy,' num2str(tbl_nme{iR,2}+1) ',' num2str(tbl_nme{iR,1}) ',LF.pre'] };
end

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_out = ejk_create_table( fcfg );

% Supp TABLE 1 - Significance (p-value)
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_pvalues.csv' ] ); end

% %%%
fcfg = [];

for iR = 1:size(tbl_nme,1)
    fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2})   ',' num2str(tbl_nme{iR,1}) ',CF.pre'] ... 
                       ['copy,' num2str(tbl_nme{iR,2})   ',' num2str(tbl_nme{iR,1}) ',LF.pre'] ...
                       ['copy,' num2str(tbl_nme{iR,2}+1) ',' num2str(tbl_nme{iR,1}) ',CF.pre'] ... 
                       ['copy,' num2str(tbl_nme{iR,2}+1) ',' num2str(tbl_nme{iR,1}) ',LF.pre'] };
end

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_pvl = ejk_create_table( fcfg );

% Supp TABLE 3 - FDR (summary)
fdr_grp_hld = cell2mat(row_lbl(:,2));
fdr_typ_hld = unique(fdr_grp_hld);

for iF = 1:numel(fdr_typ_hld)
    fdr_grp_mbm = fdr_grp_hld==fdr_typ_hld(iF);
    fdr_hld = cell2mat(tbl_pvl(fdr_grp_mbm,:));
    fdr_use{iF} = FDR(fdr_hld(:),.05);
    if isempty(fdr_use{iF}); fdr_use{iF} = 0; end
end

tbl_fdr=tbl_pvl;
for iR = 1:size(tbl_out,1)
    fdr_ind = fdr_typ_hld==row_lbl{iR,2};
    for iC = 1:size(tbl_out,2)
       if tbl_fdr{iR,iC}>fdr_use{fdr_ind}
           tbl_fdr{iR,iC}=NaN;
       end       
    end
end

% Format table
for iR = 1:size(tbl_out,1)
    for iC = 1:size(tbl_out,2)
        out_tbl{iR,iC} = num2str(roundsd(tbl_out{iR,iC},2));
        
        if tbl_pvl{iR,iC} < .01
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '**'];
        elseif tbl_pvl{iR,iC} < .05
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '*'];
        end
        
        if ~isnan(tbl_fdr{iR,iC})
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '#'];
        end
        
    end
end

out_tbl = [ col_lbl ; row_lbl(:,1) out_tbl ];

% Save Table
cell2csv( [ out_dir '/' 'AdditionalTests.csv' ], out_tbl)

%% Additional Duration Correlation Table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear out_tbl fdr_use

out_dir = [ out_put '/' 'AdditionalCorrelations' '/' ];
ejk_chk_dir(out_dir);

dta_loc = { [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/Revisions/Correlations/PreSurgical' ] ...
            [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/Revisions/Correlations/PostSurgical_Left' ] ...
            ['/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/Revisions/Correlations/PostSurgical_Right' ] };

tbl_nme = { 'Educ' 1 ; ...
            'AgeAtSurgery' 1 ; ...
            'AgeOfSeizureOnset' 1 ; ...
            'NumAEDs' 1 ; ...
            'Duration' 1 ; ...
            'Left.Hippocampus' 1 ;...
            'Right.Hippocampus' 1 ; ...
            'L.ILF' 1 ; ...
            'R.ILF' 1 ; ...
            'L.IFO' 1 ; ...
            'R.IFO' 1 ; ...
            'lh.fusiform' 1 ; ...
            'rh.fusiform' 1 };
        
col_lbl = { '' 'Pre-operative BNT' 'Pre-operative ANT' 'Post-operative BNT LTLE' 'Post-operative ANT LTLE' 'Post-operative BNT RTLE' 'Post-operative ANT RTLE' };
row_lbl = { 'Education' 1 ; ...
            'Age' 1 ; ...
            'AgeOfSeizureOnset' 1 ; ...
            'NumAEDs' 1 ; ...
            'Duration' 1 ; ...
            'Left.Hippocampus' 2 ;...
            'Right.Hippocampus' 2; ...
            'Left ILF' 2; ...
            'Right ILF' 2; ...
            'Left IFOF' 2; ...
            'Right IFOF' 2; ...
            'Left Fusiform' 2; ...
            'Right Fusiform' 2 };

% Out Table - OUTPUT (r-value)
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_rvalues.csv' ] ); end

% %%%
fcfg = [];

for iR = 1:size(tbl_nme,1)

  fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2})   ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr'] ... 
                     ['copy,' num2str(tbl_nme{iR,2})   ',' num2str(tbl_nme{iR,1}) ',ant.mem.raw.scr'] ...
                     ['copy,' num2str(tbl_nme{iR,2}+1) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr'] ... 
                     ['copy,' num2str(tbl_nme{iR,2}+1) ',' num2str(tbl_nme{iR,1}) ',ant.mem.raw.scr'] ...
                     ['copy,' num2str(tbl_nme{iR,2}+2) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr'] ... 
                     ['copy,' num2str(tbl_nme{iR,2}+2) ',' num2str(tbl_nme{iR,1}) ',ant.mem.raw.scr'] };
end

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_out = ejk_create_table( fcfg );

% Supp TABLE 1 - Significance (p-value)
clear dta_inp; for iD = 1:numel(dta_loc); dta_inp{iD} = mmil_readtext( [ dta_loc{iD} '/' 'cross_correlation_pvalues.csv' ] ); end

% %%%
fcfg = [];

for iR = 1:size(tbl_nme,1)
  fcfg.tbl(iR,:) = { ['copy,' num2str(tbl_nme{iR,2})   ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr'] ... 
                     ['copy,' num2str(tbl_nme{iR,2})   ',' num2str(tbl_nme{iR,1}) ',ant.mem.raw.scr'] ...
                     ['copy,' num2str(tbl_nme{iR,2}+1) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr'] ... 
                     ['copy,' num2str(tbl_nme{iR,2}+1) ',' num2str(tbl_nme{iR,1}) ',ant.mem.raw.scr'] ...
                     ['copy,' num2str(tbl_nme{iR,2}+2) ',' num2str(tbl_nme{iR,1}) ',bnt.raw.scr'] ... 
                     ['copy,' num2str(tbl_nme{iR,2}+2) ',' num2str(tbl_nme{iR,1}) ',ant.mem.raw.scr'] };
end

fcfg.dta = dta_inp;
fcfg.grp = grp;

tbl_pvl = ejk_create_table( fcfg );

% Supp TABLE 3 - FDR (summary)
fdr_grp_hld = cell2mat(row_lbl(:,2));
fdr_typ_hld = unique(fdr_grp_hld);

for iF = 1:numel(fdr_typ_hld)
    fdr_grp_mbm = fdr_grp_hld==fdr_typ_hld(iF);
    fdr_hld = cell2mat(tbl_pvl(fdr_grp_mbm,:));
    fdr_use{iF} = FDR(fdr_hld(:),.05);
    if isempty(fdr_use{iF}); fdr_use{iF} = 0; end
end

tbl_fdr=tbl_pvl;
for iR = 1:size(tbl_out,1)
    fdr_ind = fdr_typ_hld==row_lbl{iR,2};
    for iC = 1:size(tbl_out,2)
       if tbl_fdr{iR,iC}>fdr_use{fdr_ind}
           tbl_fdr{iR,iC}=NaN;
       end       
    end
end

% Format table
for iR = 1:size(tbl_out,1)
    for iC = 1:size(tbl_out,2)
        out_tbl{iR,iC} = num2str(roundsd(tbl_out{iR,iC},2));
        
        if tbl_pvl{iR,iC} < .01
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '**'];
        elseif tbl_pvl{iR,iC} < .05
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '*'];
        end
        
        if ~isnan(tbl_fdr{iR,iC})
            out_tbl{iR,iC} = [ out_tbl{iR,iC} '#'];
        end
        
    end
end

out_tbl = [ col_lbl ; row_lbl(:,1) out_tbl ];

% Save Table
cell2csv( [ out_dir '/' 'AdditionalDuration.csv' ], out_tbl)

%% Test Intervals & Overlap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_dir = [ out_put '/' 'Overlap' '/' ];
ejk_chk_dir(out_dir);

% Load
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Revision_Data_engel.csv'];
fcfg.dta_col = 2;
[ rev_dta, rev_dta_sbj, rev_dta_col] = ejk_dta_frm( fcfg );

% Distance
men_dst = nanmean(cell2mat(rev_dta(grp.tle_post_3T_ATLonly_all,8)));
std_dst = nanstd(cell2mat(rev_dta(grp.tle_post_3T_ATLonly_all,8)));
rng_dst(1) = min(abs(cell2mat(rev_dta(grp.tle_post_3T_ATLonly_all,8))));
rng_dst(2) = max(abs(cell2mat(rev_dta(grp.tle_post_3T_ATLonly_all,8))));

% Pre Overlap
pct_bnt_pre = sum(~isnan(cell2mat(cog_dta(grp.tle_controls_pre_3T_allSurg_all,1)))) / numel(grp.tle_controls_pre_3T_allSurg_all);
pct_ant_pre = sum(~isnan(cell2mat(cog_dta(grp.tle_controls_pre_3T_allSurg_all,2)))) / numel(grp.tle_controls_pre_3T_allSurg_all);
pct_bth_pre = sum((~isnan(cell2mat(cog_dta(grp.tle_controls_pre_3T_allSurg_all,1))) & ~isnan(cell2mat(cog_dta(grp.tle_controls_pre_3T_allSurg_all,2))))) ...
    / numel(grp.tle_controls_pre_3T_allSurg_all);

% Post Overlap
pct_bnt_pst = sum(~isnan(cell2mat(cog_dta(grp.tle_post_3T_ATLonly_all,3)))) / numel(grp.tle_post_3T_ATLonly_all);
pct_ant_pst = sum(~isnan(cell2mat(cog_dta(grp.tle_post_3T_ATLonly_all,4)))) / numel(grp.tle_post_3T_ATLonly_all);
pct_bth_pst = sum((~isnan(cell2mat(cog_dta(grp.tle_post_3T_ATLonly_all,3))) & ~isnan(cell2mat(cog_dta(grp.tle_post_3T_ATLonly_all,4))))) ...
    / numel(grp.tle_post_3T_ATLonly_all);

% Output
tbl_out = { 'Interval Mean'  men_dst*12; ...
            'Interval STD'   std_dst*12 ; ...
            'Interval Range' [num2str(rng_dst(1)*12) '-' num2str(rng_dst(2)*12)] ; ...
            'Pre BNT' pct_bnt_pre*100; ...
            'Pre ANT' pct_ant_pre*100; ...
            'Pre BNT & ANT' pct_bth_pre*100; ...
            'Post BNT' pct_bnt_pst*100; ...
            'Post ANT' pct_ant_pst*100; ...
            'Post BNT & ANT' pct_bth_pst*100 };
cell2csv([out_dir '/' 'numbers.csv'],tbl_out)
        
%% Check LTLE BNT outcomes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
out_dir = [ out_put '/' 'LTLE_Outcomes' '/' ];
ejk_chk_dir(out_dir);

% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical_wdominance.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Revision_Data_engel.csv'];
fcfg.dta_col = 2;
[ rev_dta, rev_dta_sbj, rev_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'wmparc_FA_wm_aparc_annot_QC.csv'];
fcfg.dta_col = 5;
[ wmp_dta, wmp_dta_sbj, wmp_dta_col] = ejk_dta_frm( fcfg );

% Color Scatter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grp_nme = fieldnames(grp);
iG = 8; % 8 / 

wmp_dta_one = { 'rh_fusiform' }; %
cog_dta_one = { 'bnt_raw_scr_pst'  }; % 
rev_dta_one = { 'PostInterval'      }; % 
cln_dta_one = { 'LanguageDominance' };

top_pct = 1;
cfg.col_map = { rgb('bright red') rgb('grey') rgb('bright blue') };
col_map = [];
for iC = 1:numel(cfg.col_map)-1
    col_map = [col_map ; [linspace(cfg.col_map{iC}(1),cfg.col_map{iC+1}(1),ceil(1000*top_pct/(numel(cfg.col_map)-1)))' linspace(cfg.col_map{iC}(2),cfg.col_map{iC+1}(2),ceil(1000*top_pct/(numel(cfg.col_map)-1)))' linspace(cfg.col_map{iC}(3),cfg.col_map{iC+1}(3),ceil(1000*top_pct/(numel(cfg.col_map)-1)))']; ];
end

iRO = 1;
iC = 1;

neu_dta_use = cell2mat(wmp_dta(grp.(grp_nme{iG}),  strcmpi(wmp_dta_col,wmp_dta_one{iRO})));
xdt_dta_use = cell2mat(cog_dta( grp.(grp_nme{iG}), strcmpi(cog_dta_col,cog_dta_one{iC}) ));
ydt_dta_use = cell2mat(rev_dta( grp.(grp_nme{iG}), strcmpi(rev_dta_col,rev_dta_one{iC}) ));
min_val = min(neu_dta_use);
max_val = max(neu_dta_use);
min_val = min_val - ((max_val-min_val)*.05);   

clear xdt ydt fce_col
for iCL = 1:numel(neu_dta_use)    
    xdt{iCL} = xdt_dta_use(iCL);
    ydt{iCL} = ydt_dta_use(iCL)*12;
    ndt{iCL} = neu_dta_use(iCL);
    if ~isnan(neu_dta_use(iCL))
    fce_col{iCL} = col_map(ceil(((neu_dta_use(iCL)-min_val) / (max_val - min_val))*1000)-1,:);
    else
        fce_col{iCL} = rgb('white');
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.xdt     = xdt([1:7 9:end]);
fcfg.ydt     = num2cell(cellfun(@abs,ydt([1:7 9:end])));

fcfg.fce_col = fce_col;
fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));

fcfg.xlb = { cog_dta_one{iC} };
fcfg.ylb = { [rev_dta_one{iC} '_' 'months']  };

fcfg.hln     = -1.28;
fcfg.hln_col = rgb('red');

fcfg.vln     = -1.28;
fcfg.vln_col = rgb('red');

fcfg.jtr = 0;
fcfg.xlm = [ -5 5 ];
fcfg.ylm = [ 0 24 ];

fcfg.out_dir = out_dir;
fcfg.out_nme = [ grp_nme{iG} '_' cog_dta_one{iC} '_BY_' rev_dta_one{iC} '_' wmp_dta_one{iRO} '_no_out_24' ];

ejk_scatter(fcfg)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.xdt     = ndt([1:7 9:end]);
fcfg.ydt     = num2cell(cellfun(@abs,ydt([1:7 9:end])));

fcfg.fce_col = fce_col;
fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));

fcfg.xlb = { wmp_dta_one{iRO}};
fcfg.ylb = { [rev_dta_one{iC} '_' 'months']  };

fcfg.hln     = -1.28;
fcfg.hln_col = rgb('red');

fcfg.vln     = -1.28;
fcfg.vln_col = rgb('red');

fcfg.jtr = 0;
% fcfg.xlm = [ -5 5 ];
fcfg.ylm = [ 0 24 ];

fcfg.out_dir = out_dir;
fcfg.out_nme = [ grp_nme{iG} '_' wmp_dta_one{iRO} '_BY_' rev_dta_one{iC} '_no_out' ];

ejk_scatter(fcfg)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.xdt     = xdt([1:7 9:end]);
fcfg.ydt     = ndt([1:7 9:end]);

fcfg.fce_col = fce_col;
fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));

fcfg.xlb = { cog_dta_one{iC}};
fcfg.ylb = { wmp_dta_one{iRO}  };

fcfg.hln     = -1.28;
fcfg.hln_col = rgb('red');

fcfg.vln     = -1.28;
fcfg.vln_col = rgb('red');

fcfg.jtr = 0;
fcfg.xlm = [ -5 5 ];
% fcfg.ylm = [ 0 24 ];

fcfg.out_dir = out_dir;
fcfg.out_nme = [ grp_nme{iG} '_' wmp_dta_one{iRO} '_BY_' cog_dta_one{iC} '_no_out' ];

ejk_scatter(fcfg)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.xdt     = xdt;
fcfg.ydt     = num2cell(cellfun(@abs,ydt));

fcfg.fce_col = fce_col;
fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));

fcfg.xlb = { cog_dta_one{iC} };
fcfg.ylb = { [rev_dta_one{iC} '_' 'months']  };

fcfg.hln     = -1.28;
fcfg.hln_col = rgb('red');

fcfg.vln     = -1.28;
fcfg.vln_col = rgb('red');

fcfg.jtr = 0;
fcfg.xlm = [ -5 5 ];
% fcfg.ylm = [ 0 24 ];

fcfg.out_dir = out_dir;
fcfg.out_nme = [ grp_nme{iG} '_' cog_dta_one{iC} '_BY_' rev_dta_one{iC} '_' wmp_dta_one{iRO} '_24' ];

ejk_scatter(fcfg)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.xdt     = ndt;
fcfg.ydt     = num2cell(cellfun(@abs,ydt));

fcfg.fce_col = fce_col;
fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));

fcfg.xlb = { wmp_dta_one{iRO}};
fcfg.ylb = { [rev_dta_one{iC} '_' 'months']  };

fcfg.hln     = -1.28;
fcfg.hln_col = rgb('red');

fcfg.vln     = -1.28;
fcfg.vln_col = rgb('red');

fcfg.jtr = 0;
% fcfg.xlm = [ -5 5 ];
% fcfg.ylm = [ 0 24 ];

fcfg.out_dir = out_dir;
fcfg.out_nme = [ grp_nme{iG} '_' wmp_dta_one{iRO} '_BY_' rev_dta_one{iC} ];

ejk_scatter(fcfg)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.xdt     = xdt;
fcfg.ydt     = ndt;

fcfg.fce_col = fce_col;
fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));

fcfg.xlb = { cog_dta_one{iC}};
fcfg.ylb = { wmp_dta_one{iRO}  };

fcfg.hln     = -1.28;
fcfg.hln_col = rgb('red');

fcfg.vln     = -1.28;
fcfg.vln_col = rgb('red');

fcfg.jtr = 0;
fcfg.xlm = [ -5 5 ];
% fcfg.ylm = [ 0 24 ];

fcfg.out_dir = out_dir;
fcfg.out_nme = [ grp_nme{iG} '_' wmp_dta_one{iRO} '_BY_' cog_dta_one{iC} ];

ejk_scatter(fcfg)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ rvl, pvl] = corrcoef( cell2mat(xdt([1:7 9:end])), cellfun(@abs,ydt([1:7 9:end])), 'Rows', 'pairwise' );
[ rvl, pvl] = corr( cell2mat(xdt([1:7 9:end]))', cellfun(@abs,ydt([1:7 9:end]))', 'Type', 'Spearman', 'Rows', 'pairwise' );

[ rvl_neu_cog, pvl_neu_cog] = corr(        neu_dta_use, xdt_dta_use, 'Type', 'Spearman', 'Rows', 'pairwise' );
[ rvl_cog_int, pvl_cog_int] = corr(        ydt_dta_use, xdt_dta_use, 'Type', 'Spearman', 'Rows', 'pairwise' );
[ rvl_neu_prt, pvl_neu_prt] = partialcorr( neu_dta_use, xdt_dta_use, ydt_dta_use, 'Type', 'Spearman', 'Rows', 'pairwise' );
out_csv(1,:) = num2cell([ pvl_neu_cog pvl_cog_int pvl_neu_prt ]);

[ rvl_neu_cog, pvl_neu_cog] = corr(        neu_dta_use([1:7 9:end]), xdt_dta_use([1:7 9:end]), 'Type', 'Spearman', 'Rows', 'pairwise' );
[ rvl_cog_int, pvl_cog_int] = corr(        ydt_dta_use([1:7 9:end]), xdt_dta_use([1:7 9:end]), 'Type', 'Spearman', 'Rows', 'pairwise' );
[ rvl_neu_prt, pvl_neu_prt] = partialcorr( neu_dta_use([1:7 9:end]), xdt_dta_use([1:7 9:end]), ydt_dta_use([1:7 9:end]), 'Type', 'Spearman', 'Rows', 'pairwise' );
out_csv(2,:) = num2cell([ pvl_neu_cog pvl_cog_int pvl_neu_prt ]);

cell2csv([ out_dir '/' 'pvals.csv'], out_csv);

[ rvl_cog_int, pvl_cog_int] = corr(        ydt_dta_use([1:7 9:end]), xdt_dta_use([1:7 9:end]), 'Type', 'Spearman', 'Rows', 'pairwise' );
[ rvl_neu_prt, pvl_neu_prt] = partialcorr( ydt_dta_use([1:7 9:end]), xdt_dta_use([1:7 9:end]), neu_dta_use([1:7 9:end]), 'Type', 'Spearman', 'Rows', 'pairwise' );
[pvl_cog_int pvl_neu_prt]

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iRO = 1;
iC = 1;

neu_dta_use = cell2mat(wmp_dta(grp.(grp_nme{iG}),  strcmpi(wmp_dta_col,wmp_dta_one{iRO})));
xdt_dta_use = cell2mat(cog_dta( grp.(grp_nme{iG}), strcmpi(cog_dta_col,cog_dta_one{iC}) ));
cln_dta_use = cln_dta( grp.(grp_nme{iG}), strcmpi(cln_dta_col,cln_dta_one{iC}) );

clear xdt ydt fce_col
for iCL = 1:numel(neu_dta_use)    
    xdt{iCL} = xdt_dta_use(iCL);
    ydt{iCL} = neu_dta_use(iCL);
    if strcmpi(cln_dta_use(iCL),'Atypical')
        fce_col{iCL} = rgb('mauve');
    elseif strcmpi(cln_dta_use(iCL),'Typical')
        fce_col{iCL} = rgb('blue');
    else
        fce_col{iCL} = rgb('white');
    end
end

fcfg = [];

fcfg.xdt     = xdt;
fcfg.ydt     = ydt;

fcfg.fce_col = fce_col;
fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));

fcfg.xlb = { cog_dta_one{iC} };
fcfg.ylb = { wmp_dta_one{iRO}  };

fcfg.vln     = -1.28;
fcfg.vln_col = rgb('red');

fcfg.jtr = 0;
fcfg.xlm = [ -5 5 ];

fcfg.out_dir = out_dir;
fcfg.out_nme = [ grp_nme{iG} '_' cog_dta_one{iC} '_BY_' wmp_dta_one{iRO} '_typical_atypical' ];

ejk_scatter(fcfg)

%% Entorhinal Figure
fcfg = [];

fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc'; 
fcfg.fsr_nme = 'fsaverage';

fcfg.roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer';
fcfg.prc_nme = '.aparc.annot';

fcfg.inc_reg = { 'entorhinal' };

fcfg.sph = { 'lh' 'rh' };
fcfg.sph_vew = { 'lat' 'ven' 'med' };

fcfg.out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/atl_nme/Submission/NeuroimageClinical/Revision/Additional/NewFig';
fcfg.out_nme  = 'entorhinal_roi';

ejk_roi_plot(fcfg);


